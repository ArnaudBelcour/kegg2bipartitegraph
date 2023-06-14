# Copyright (C) 2021-2023 Arnaud Belcour - Inria Dyliss
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

from bioservices import version as bioservices_version
from bioservices import KEGG, UniProt

from kegg2bipartitegraph.utils import is_valid_dir, get_rest_uniprot_release
from kegg2bipartitegraph.reference import get_kegg_database_version
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

# Create KEGG instance of bioservices.KEEG.
KEGG_BIOSERVICES = KEGG()
UNIPROT_BIOSERVICES = UniProt(verbose=False)

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
    print(data)
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
    print(len(protein_to_search_on_uniprots))
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

def old_create_organism_sbml_cobrapy(base_filename, kegg_sbml_model_path, taxon_reactions, clust_sbml_output_folder_path):
    """ Unused function only kept at this moment for archiving.
    """
    from cobra.io import read_sbml_model, write_sbml_model
    from cobra import Model
    kegg_model = read_sbml_model(kegg_sbml_model_path)
    kegg_model = read_sbml_model(kegg_sbml_model_path)

    species_model = Model(base_filename)

    # Map KEGG reaction from KEGG SBML model to taxon SBML.
    sbml_reactions = []
    for reaction in kegg_model.reactions:
        reaction_id = reaction.id.replace('R_','')
        if reaction_id in taxon_reactions:
            reaction.gene_reaction_rule = '( ' + ' or '.join(taxon_reactions[reaction_id]) + ' )'
            sbml_reactions.append(reaction)

    species_model.add_reactions(sbml_reactions)

    # Create file if there is at least 1 reaction.
    if len(species_model.reactions) > 0:
        # Create SBML file.
        sbml_output_file_path = os.path.join(clust_sbml_output_folder_path, base_filename+'.sbml')
        write_sbml_model(species_model, sbml_output_file_path)
    else:
        logger.info('|EsMeCaTa|kegg| No reactions in model for {0}, no SBML file will be created.'.format(base_filename))


def create_organism_sbml(kegg_sbml_model_path, taxon_reactions):
    """Create a SBML model from the KEGG reference model and reference found for the taxon.

    Args:
        kegg_sbml_model_path (str): KEGG code for an organism
        taxon_reactions (dict): reaction ID as key and associated genes as value
    """
    # Read the reference KEGG sbml file.
    # Use it to create the organism sbml file.
    reader = libsbml.SBMLReader()
    kegg_document = reader.readSBML(kegg_sbml_model_path)
    kegg_model = kegg_document.getModel()

    # Remove all the gene products from the model.
    model_fbc = kegg_model.getPlugin('fbc')
    model_fbc.setStrict(True)

    remove_gene_products = []
    for gene_product in model_fbc.getListOfGeneProducts():
        remove_gene_products.append(gene_product.id)
    for gene_product in remove_gene_products:
        model_fbc.removeGeneProduct(gene_product)

    # Keep only the reactions find during the reconstruction process.
    genes = []
    remove_reactions = []
    kept_metabolites = []
    for reaction in kegg_model.getListOfReactions():
        if reaction.id in taxon_reactions:
            # Add the gene associated with the organism.
            r_fbc: "libsbml.FbcReactionPlugin" = reaction.getPlugin("fbc")
            gpr_association = r_fbc.createGeneProductAssociation()

            for gene in taxon_reactions[reaction.id]:
                if gene not in genes:
                    gene_prod = model_fbc.createGeneProduct()
                    gene_prod.setId(gene), 'add gene %s' %gene
                    gene_prod.setName(gene)
                    gene_prod.setLabel(gene)
                    genes.append(gene)
            gpr = ' or '.join(taxon_reactions[reaction.id])
            gpr_association.setAssociation(gpr, True, True)
            kept_metabolites.extend([i.species for i in reaction.getListOfReactants()])
            kept_metabolites.extend([i.species for i in reaction.getListOfProducts()])
        else:
            remove_reactions.append(reaction.id)

    # Remove reactions not found in organism.
    for reaction_id in remove_reactions:
        kegg_model.removeReaction(reaction_id)
    # Remove metabolites not found in organism.
    remove_metabolites = set([m.id for m in kegg_model.getListOfSpecies()]) - set(kept_metabolites)
    for metabolite_id in remove_metabolites:
        kegg_model.removeSpecies(metabolite_id)

    return kegg_document, kegg_model


def create_draft_networks(input_folder, output_folder, mapping_ko=False, recreate_kegg=None, remove_ubiquitous=True):
    """From the output folder of 'esmecata annotation' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        input_folder (str): path to the output folder of esmecata annotation
        output_folder (str): path to the output folder
        mapping_ko (bool): option to use KO ID to retrieve reactions
        recreate_kegg (bool): option to recreate KEGG model by guerying KEGG server
    """
    starttime = time.time()
    logger.info('|EsMeCaTa|kegg| Begin KEGG metabolism mapping.')

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder, 'mapping_ko': mapping_ko}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__

    if mapping_ko:
        kegg2bipartitegraph_esmecata_metadata = get_rest_uniprot_release(options)
        kegg2bipartitegraph_esmecata_metadata['kegg_release_number'] = get_kegg_database_version()
    else:
        kegg2bipartitegraph_esmecata_metadata = {}
        kegg2bipartitegraph_esmecata_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    data_kegg_model_path = DATA_ROOT
    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(data_kegg_model_path, 'kegg_model')
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_esmecata_metadata['reference_path'] = data_kegg_model_path

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')

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

    # Retrieve proteins and annotations from esmecata annotation folder.
    annotation_reference_folder_path = os.path.join(input_folder, 'annotation_reference')

    if mapping_ko is True:
        clust_ko_output_folder_path = ko_output_folder_path

    clust_pathways_output_folder_path = pathways_output_folder_path
    is_valid_dir(clust_pathways_output_folder_path)

    clust_sbml_output_folder_path = sbml_output_folder_path
    is_valid_dir(clust_sbml_output_folder_path)

    for annot_file in os.listdir(annotation_reference_folder_path):
        annot_file_path = os.path.join(annotation_reference_folder_path, annot_file)

        base_file = os.path.basename(annot_file_path)
        base_filename = os.path.splitext(base_file)[0]

        # Extract protein IDs and EC number from anntotation reference folder.
        protein_ec_numbers = {}
        protein_clusters = {}
        with open(annot_file_path, 'r') as open_annot_file_path:
            csvreader = csv.reader(open_annot_file_path, delimiter='\t')
            next(csvreader)
            for line in csvreader:
                protein_id = line[0]
                protein_cluster = line[1].split(',')
                ec_numbers = line[4].split(',')
                protein_ec_numbers[protein_id] = ec_numbers
                protein_clusters[protein_id] = protein_cluster

        taxon_reactions = {}
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
            logger.info('|EsMeCaTa|kegg| Added {0} reactions from {1} KO for taxon {2}.'.format(len(set(ko_added_reactions)), len(set(all_kos)), base_filename))

        # Use EC found to be associated to reference protein to retrieve KEEG reaction.
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
        logger.info('|EsMeCaTa|kegg| Added {0} reactions from {1} EC for taxon {2}.'.format(len(set(ec_added_reactions)), len(set(all_ecs)), base_filename))

        total_added_reactions = list(taxon_reactions.keys())
        if mapping_ko is True:
            logger.info('|EsMeCaTa|kegg| A total of {0} unique reactions are added from EC and KO for taxon {1}.'.format(len(total_added_reactions), base_filename))

        # Create pathway file contening pathway with reactions in the taxon.
        pathways_output_file_path = os.path.join(clust_pathways_output_folder_path, base_filename+'.tsv')
        with open(pathways_output_file_path, 'w') as open_pathways_output_file_path:
            csvwriter = csv.writer(open_pathways_output_file_path, delimiter='\t')
            csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_completion_ratio', 'pathway_reaction_in_taxon', 'pathway_reaction'])
            for pathway in kegg_pathways:
                pathway_reactions = kegg_pathways[pathway][1]
                pathway_reaction_in_taxon = set(pathway_reactions).intersection(set(total_added_reactions))
                if len(pathway_reaction_in_taxon) > 0:
                    pathway_name = kegg_pathways[pathway][0]
                    pathway_completion_ratio = len(pathway_reaction_in_taxon) / len(pathway_reactions)
                    csvwriter.writerow([pathway, pathway_name, pathway_completion_ratio, ','.join(pathway_reaction_in_taxon), ','.join(pathway_reactions)])

        # Create module file contening module with reactions in the taxon.
        modules_output_file_path = os.path.join(modules_output_folder_path, base_filename+'.tsv')
        with open(modules_output_file_path, 'w') as open_modules_output_file_path:
            csvwriter = csv.writer(open_modules_output_file_path, delimiter='\t')
            csvwriter.writerow(['module_id', 'module_name', 'module_completion_ratio', 'module_reaction_in_taxon', 'module_reaction'])
            for module in kegg_modules:
                module_reactions = kegg_modules[module][1]
                module_reaction_in_taxon = set(module_reactions).intersection(set(total_added_reactions))
                if len(module_reaction_in_taxon) > 0:
                    module_name = kegg_modules[module][0]
                    module_completion_ratio = len(module_reaction_in_taxon) / len(module_reactions)
                    csvwriter.writerow([module, module_name, module_completion_ratio, ','.join(module_reaction_in_taxon), ','.join(module_reactions)])

        kegg_document, kegg_model = create_organism_sbml(kegg_sbml_model_path, taxon_reactions)

        # Create file if there is at least 1 reaction.
        if len(kegg_model.getListOfReactions()) > 0:
            # Create SBML file.
            sbml_output_file_path = os.path.join(sbml_output_folder_path, base_filename+'.sbml')
            libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
            logger.info('|kegg2bipartitegraph|organism| Network of {0} contains {1} reactions and {2} metabolites.'.format(base_filename, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
        else:
            logger.info('|kegg2bipartitegraph|organism| No reactions in model for {0}, no SBML file will be created.'.format(base_filename))

    clust_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    compute_stat_kegg(clust_sbml_output_folder_path, clust_stat_file)

    endtime = time.time()

    duration = endtime - starttime
    kegg2bipartitegraph_esmecata_metadata['kegg2bipartitegraph_esmecata_duration'] = duration
    kegg2bipartitgraph_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_esmecata_kegg.json')
    with open(kegg2bipartitgraph_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_esmecata_metadata, ouput_file, indent=4)
    logger.info('|EsMeCaTa|kegg| Draft networks creation complete.')
