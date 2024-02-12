# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import csv
import json
import logging
import os
import libsbml
import time
import sys
import urllib.parse
import urllib.request

from networkx import __version__ as networkx_version
from bioservices import version as bioservices_version
from bioservices import KEGG, UniProt

from kegg2bipartitegraph.utils import is_valid_dir, get_rest_uniprot_release, write_pathway_file, write_module_file
from kegg2bipartitegraph.reference import get_kegg_database_version
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg, get_go_to_ec
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def chunks(elements, n):
    """Yield successive n-sized chunks from list.
    Form: https://stackoverflow.com/a/312464
    Args:
        elements (list): list of elements (proteins or proteomes) to be split in chunks of size n
        n (int): size of the chunks
    Returns:
        list: list of elements of size n
    """
    for i in range(0, len(elements), n):
        yield elements[i:i + n]


def query_uniprot_bioservices(protein_queries):
    """REST query to get annotation from proteins.
    Args:
        protein_queries (list): list of proteins.
    Returns:
        data (dict): dictionary returned by bioservices containing annotation
    """
    data = UNIPROT_BIOSERVICES.mapping(fr='UniProtKB_AC-ID', to='KEGG', query=protein_queries,
                                        progress=True)
    return data


def query_uniprot_kegg_rest(protein_to_search_on_uniprots, output_dict):
    """Query UniProt with REST to find KEGG gene ID

    Args:
        protein_to_search_on_uniprots (set): set of proteins not already annotated that will need UniProt queries
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value

    Returns:
        output_dict (dict): annotation dict: protein as key and KEGG gene ID as value
    """
    # The limit of 15 000 proteins per query comes from the help of Uniprot (inferior to 20 000):
    # https://www.uniprot.org/help/uploadlists
    if len(protein_to_search_on_uniprots) < 5000:
        protein_queries = ','.join(protein_to_search_on_uniprots)
        tmp_output_dict = query_uniprot_bioservices(protein_queries)
        output_dict.update(tmp_output_dict)
        time.sleep(1)
    else:
        protein_chunks = chunks(list(protein_to_search_on_uniprots), 5000)
        for chunk in protein_chunks:
            protein_queries = ','.join(chunk)
            tmp_output_dict = query_uniprot_bioservices(protein_queries)
            output_dict.update(tmp_output_dict)
            time.sleep(1)

    return output_dict


def map_protein_to_KO_id(proteins_ids):
    """From a list of UniProt protein IDs, retrieved the corresponding KEGG Orthologs.

    Args:
        proteins_ids (list): List of Uniprot protein IDs

    Returns:
        protein_ko_mapping (dict): mapping between UniProt protein ID as key and KO ID as value
    """
    protein_mapping = {}
    protein_mapping = query_uniprot_kegg_rest(proteins_ids, protein_mapping)
    kegg_genes = [gene for prot_id in protein_mapping for gene in protein_mapping[prot_id]]

    chunk_kegg_genes = chunks(kegg_genes, 99)

    mapping_gene_kegg_ko = {}
    for chunk_kegg_gene in chunk_kegg_genes:
        str_chunk_kegg_gene = '+'.join(chunk_kegg_gene)
        response_text = KEGG_BIOSERVICES.link('ko', str_chunk_kegg_gene)
        if response_text != '':
            csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
        else:
            csvreader = []
        for line in csvreader:
            kegg_gene = line[0]
            ko = line[1]
            mapping_gene_kegg_ko[kegg_gene] = ko

    protein_ko_mapping = {}
    for protein in protein_mapping:
        kegg_genes = protein_mapping[protein]
        for kegg_gene in kegg_genes:
            if kegg_gene in mapping_gene_kegg_ko:
                ko = mapping_gene_kegg_ko[kegg_gene]
                protein_ko_mapping[protein] = ko

    return protein_ko_mapping


def create_esmecata_network(input_folder, output_folder, mapping_ko=False, reference_folder=False):
    """From the output folder of 'esmecata annotation' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        input_folder (str): path to the output folder of esmecata annotation
        output_folder (str): path to the output folder
        mapping_ko (bool): option to use KO ID to retrieve reactions
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|esmecata| Begin KEGG metabolism mapping.')

    if reference_folder is not False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|esmecata| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|esmecata| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder, 'mapping_ko': mapping_ko}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__
    options['tool_dependencies']['python_package']['networkx'] = networkx_version

    if mapping_ko:
        # Create KEGG instance of bioservices.KEEG in this function to avoid trying to connect to KEGG with offline mode.
        global KEGG_BIOSERVICES
        KEGG_BIOSERVICES = KEGG()
        global UNIPROT_BIOSERVICES
        UNIPROT_BIOSERVICES = UniProt(verbose=False)
        kegg2bipartitegraph_esmecata_metadata = get_rest_uniprot_release(options)
        kegg2bipartitegraph_esmecata_metadata['kegg_release_number'] = get_kegg_database_version()
    else:
        kegg2bipartitegraph_esmecata_metadata = {}
        kegg2bipartitegraph_esmecata_metadata['tool_options'] = options

    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_esmecata_metadata['reference_path'] = kegg_model_path

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_json_model_path = os.path.join(kegg_model_path, 'kegg_metadata.json')
    kegg_json_hierarchy_path = os.path.join(kegg_model_path, 'kegg_hierarchy.json')
    ec_to_gos_path = os.path.join(kegg_model_path, 'ec_to_gos.tsv')

    # Read the reference KEGG sbml file.
    # Use it to create the organism sbml file.
    reader = libsbml.SBMLReader()
    reference_kegg_document = reader.readSBML(kegg_sbml_model_path)
    reference_kegg_model = reference_kegg_document.getModel()
    reference_model_fbc = reference_kegg_model.getPlugin('fbc')

    reference_reactions = {reaction.id: reaction for reaction in reference_kegg_model.getListOfReactions()}
    reference_species = reference_kegg_model.getListOfSpecies()
    reference_groups = reference_kegg_model.getPlugin("groups")

    with open(kegg_json_model_path, 'r') as input_metadata_json:
        json_data = json.load(input_metadata_json)

    kegg2bipartitegraph_esmecata_metadata['json_reference_metadata'] = json_data

    kegg_pathways = {}
    with open(kegg_pathways_path, 'r') as open_kegg_pathways_path:
        csvreader = csv.reader(open_kegg_pathways_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            pathway_id = line[0]
            pathway_name = line[1]
            pathway_reactions = line[2].split(',')
            kegg_pathways[pathway_id] = (pathway_name, pathway_reactions)

    kegg_modules = {}
    with open(kegg_modules_path, 'r') as open_kegg_modules_path:
        csvreader = csv.reader(open_kegg_modules_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            moudle_id = line[0]
            module_name = line[1]
            module_reactions = line[3].split(',')
            kegg_modules[moudle_id] = (module_name, module_reactions)

    # Create SBML output folder.
    sbml_output_folder_path = os.path.join(output_folder, 'sbml')
    is_valid_dir(sbml_output_folder_path)

    # Create graphml output folder.
    graphml_output_folder_path = os.path.join(output_folder, 'graphml')
    is_valid_dir(graphml_output_folder_path)

    # Create KO output folder.
    ko_output_folder_path = os.path.join(output_folder, 'ko')
    is_valid_dir(ko_output_folder_path)

    # Create pathways annotated output folder.
    pathways_output_folder_path = os.path.join(output_folder, 'pathways')
    is_valid_dir(pathways_output_folder_path)

    # Create modules annotated output folder.
    modules_output_folder_path = os.path.join(output_folder, 'modules')
    is_valid_dir(modules_output_folder_path)

    ko_to_reactions, ec_to_reactions = retrieve_mapping_dictonaries(kegg_rxn_mapping_path)
    go_to_ecs = get_go_to_ec(ec_to_gos_path)

    # Retrieve proteins and annotations from esmecata annotation folder.
    annotation_reference_folder_path = os.path.join(input_folder, 'annotation_reference')

    if mapping_ko is True:
        clust_ko_output_folder_path = ko_output_folder_path

    clust_pathways_output_folder_path = pathways_output_folder_path
    is_valid_dir(clust_pathways_output_folder_path)

    clust_sbml_output_folder_path = sbml_output_folder_path
    is_valid_dir(clust_sbml_output_folder_path)

    # Extract class modules
    with open(kegg_json_hierarchy_path, 'r') as open_hierarchy_json:
        hierarchy_json_data = json.load(open_hierarchy_json)
    hierarchy_json_data_module = hierarchy_json_data['module']['Pathway modules']

    class_representation = {}
    for module_class in hierarchy_json_data_module:
        for module_subclass in hierarchy_json_data_module[module_class]:
            class_representation[module_subclass] = set(hierarchy_json_data_module[module_class][module_subclass])

    class_modules = [class_module for class_module in class_representation]
    class_data = []

    stat_metabolic_networks = {}
    for annot_file in os.listdir(annotation_reference_folder_path):
        annot_file_path = os.path.join(annotation_reference_folder_path, annot_file)

        base_file = os.path.basename(annot_file_path)
        base_filename = os.path.splitext(base_file)[0]

        taxon_reactions = {}
        # Extract protein IDs and EC number from anntotation reference folder.
        protein_ec_numbers = {}
        protein_go_numbers = {}
        protein_kegg_reactions = {}
        protein_clusters = {}
        with open(annot_file_path, 'r') as open_annot_file_path:
            csvreader = csv.DictReader(open_annot_file_path, delimiter='\t')
            for line in csvreader:
                protein_id = line['protein_cluster']
                protein_cluster = line['cluster_members'].split(',')
                ec_numbers = line['EC'].split(',')
                go_numbers = line['GO'].split(',')
                if 'KEGG_reaction' in line:
                    kegg_reactions = line['KEGG_reaction'].split(',')
                else:
                    kegg_reactions = []

                protein_ec_numbers[protein_id] = ec_numbers
                protein_go_numbers[protein_id] = [go for go in go_numbers if go != '']
                protein_kegg_reactions[protein_id] = kegg_reactions
                protein_clusters[protein_id] = protein_cluster

        # If mapping KO option is used, search for KO terms associated to proteins.
        if mapping_ko is True:
            protein_to_maps = set([protein_id for protein_cluster in protein_clusters for protein_id in protein_clusters[protein_cluster]])
            protein_ko_mapping = map_protein_to_KO_id(protein_to_maps)

            ko_output_file_path = os.path.join(clust_ko_output_folder_path, base_filename+'.tsv')
            with open(ko_output_file_path, 'w') as output_file:
                csvwriter = csv.writer(output_file, delimiter='\t')
                csvwriter.writerow(['protein', 'KO'])
                for protein_cluster in protein_clusters:
                    for protein_id in protein_clusters[protein_cluster]:
                        if protein_id in protein_ko_mapping:
                            ko = protein_ko_mapping[protein_id]
                        else:
                            ko = ''
                        csvwriter.writerow([protein_id, ko])

            # Keep KO and propagate their reaction only if they appear in all protein of the cluster.
            ko_added_reactions = []
            all_kos = []
            for protein_cluster in protein_clusters:
                ko_ids = [[protein_ko_mapping[protein_id]] if protein_id in protein_ko_mapping else [] for protein_id in protein_clusters[protein_cluster]]
                protein_all_ko_ids = set([ko_id for subko_ids in ko_ids for ko_id in subko_ids])
                all_kos.extend(list(protein_all_ko_ids))
                keep_kos = [ko_id for ko_id in protein_all_ko_ids if sum(subko_ids.count(ko_id) for subko_ids in ko_ids) >= 1 * len(protein_clusters[protein_cluster])]
                for ko_id in keep_kos:
                    ko_id = ko_id.replace('ko:', '')
                    if ko_id in ko_to_reactions:
                        reaction_ids = ko_to_reactions[ko_id]
                        for reaction_id in reaction_ids:
                            ko_added_reactions.append(reaction_id)
                            if reaction_id not in taxon_reactions:
                                taxon_reactions[reaction_id] = [protein_cluster]
                            else:
                                taxon_reactions[reaction_id].append(protein_cluster)
            logger.info('|kegg2bipartitegraph|esmecata| Added {0} reactions from {1} KO for taxon {2}.'.format(len(set(ko_added_reactions)), len(set(all_kos)), base_filename))

        # Extract KEGG reacitons form results.
        kegg_reaction_added_reactions = []
        all_kegg_reactions = []
        for protein in protein_kegg_reactions:
            kegg_reactions = protein_kegg_reactions[protein]
            all_kegg_reactions.extend(kegg_reactions)
            for kegg_reaction in kegg_reactions:
                kegg_reaction_added_reactions.append(kegg_reaction)
                if kegg_reaction not in taxon_reactions:
                    taxon_reactions[kegg_reaction] = [protein]
                else:
                    taxon_reactions[kegg_reaction].append(protein)
        logger.info('|kegg2bipartitegraph|esmecata| Added {0} reactions from {1} KEGG reactions for taxon {2}.'.format(len(set(kegg_reaction_added_reactions)), len(set(all_kegg_reactions)), base_filename))

        # Extract GO terms and map them to KEGG reaction.
        go_added_reactions = []
        all_gos = []
        for protein in protein_go_numbers:
            go_ids = protein_go_numbers[protein]
            all_gos.extend(go_ids)
            for go_id in go_ids:
                if go_id in go_to_ecs:
                    ecs = go_to_ecs[go_id]
                    for ec in ecs:
                        if ec in ec_to_reactions:
                            reaction_ids = ec_to_reactions[ec]
                            for reaction_id in reaction_ids:
                                go_added_reactions.append(reaction_id)
                                if reaction_id not in taxon_reactions:
                                    taxon_reactions[reaction_id] = [protein]
                                else:
                                    taxon_reactions[reaction_id].append(protein)
        logger.info('|kegg2bipartitegraph|esmecata| Added {0} reactions from {1} GO for taxon {2}.'.format(len(set(go_added_reactions)), len(set(all_gos)), base_filename))

        # Extract EC numbers and map them to KEGG reaction.
        ec_added_reactions = []
        all_ecs = []
        for protein in protein_ec_numbers:
            ec_ids = protein_ec_numbers[protein]
            all_ecs.extend(ec_ids)
            for ec_id in ec_ids:
                if ec_id in ec_to_reactions:
                    reaction_ids = ec_to_reactions[ec_id]
                    for reaction_id in reaction_ids:
                        ec_added_reactions.append(reaction_id)
                        if reaction_id not in taxon_reactions:
                            taxon_reactions[reaction_id] = [protein]
                        else:
                            taxon_reactions[reaction_id].append(protein)
        logger.info('|kegg2bipartitegraph|esmecata| Added {0} reactions from {1} EC for taxon {2}.'.format(len(set(ec_added_reactions)), len(set(all_ecs)), base_filename))

        total_added_reactions = list(taxon_reactions.keys())
        if mapping_ko is True:
            logger.info('|kegg2bipartitegraph|esmecata| A total of {0} unique reactions are added from EC and KO for taxon {1}.'.format(len(total_added_reactions), base_filename))

        # Create pathway file contening pathway with reactions in the taxon.
        pathways_output_file_path = os.path.join(clust_pathways_output_folder_path, base_filename+'.tsv')
        organism_pathways = write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions)

        # Create module file contening module with reactions in the taxon.
        modules_output_file_path = os.path.join(modules_output_folder_path, base_filename+'.tsv')
        organism_modules = write_module_file(kegg_modules, modules_output_file_path, total_added_reactions)
        class_data.append([base_filename, *[len(class_representation[class_module].intersection(set(organism_modules))) for class_module in class_modules]])

        kegg_document, kegg_model = create_sbml_from_kegg_reactions(base_filename, reference_reactions, reference_species, reference_groups, taxon_reactions, organism_pathways, organism_modules)

        stat_metabolic_networks[base_filename] = (len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies()))

        # Create file if there is at least 1 reaction.
        if len(kegg_model.getListOfReactions()) > 0:
            # Create SBML file.
            sbml_output_file_path = os.path.join(sbml_output_folder_path, base_filename+'.sbml')
            libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
            graphml_output_file_path = os.path.join(graphml_output_folder_path, base_filename+'.graphml')
            sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
            logger.info('|kegg2bipartitegraph|esmecata| Network of {0} contains {1} reactions and {2} metabolites.'.format(base_filename, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
        else:
            logger.info('|kegg2bipartitegraph|esmecata| No reactions in model for {0}, no SBML file will be created.'.format(base_filename))

    module_class_file = os.path.join(output_folder, 'module_class.tsv')
    with open(module_class_file, 'w') as open_module_class_file:
        csvwriter = csv.writer(open_module_class_file, delimiter='\t')
        csvwriter.writerow(['observation_name', *class_modules])
        for class_dat in class_data:
            csvwriter.writerow(class_dat)

    kegg_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    with open(kegg_stat_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_reactions', 'Number_metabolites'])
        for observation_name in stat_metabolic_networks:
            csvwriter.writerow([observation_name, stat_metabolic_networks[observation_name][0], stat_metabolic_networks[observation_name][1]])

    endtime = time.time()

    duration = endtime - starttime
    kegg2bipartitegraph_esmecata_metadata['kegg2bipartitegraph_esmecata_duration'] = duration
    kegg2bipartitgraph_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_esmecata_kegg.json')
    with open(kegg2bipartitgraph_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_esmecata_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|esmecata| Draft networks creation complete.')
