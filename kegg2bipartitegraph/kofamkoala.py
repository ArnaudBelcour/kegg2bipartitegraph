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
import logging
import os
import json
import time
import pandas as pd
import re
import sys
import libsbml
import urllib.request

from networkx import __version__ as networkx_version
from bioservices import version as bioservices_version
from kegg2bipartitegraph.utils import is_valid_dir, write_pathway_file, write_module_file
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def read_kofam_koala_txt(kofam_koala_result_file):
    """Read a kofam koala annotation file and retrieve KO with gene ID.

    Args:
        kofam_koala_result_file (str): path to kofam koala annotation file
    Returns:
        dict: dict of genes and their annotations as {gene1:[KO1, KO2]}
    """
    gene_to_kos = {}

    ko_regex = r'K\d{5}'
    with open(kofam_koala_result_file) as open_kofam_file:
        for line in open_kofam_file.readlines():
            if '*' in line:
                # * indicates result above threshold and is followed by gene ID.
                gene_id = line.split('* ')[1].split(' ')[0]
                if re.search(ko_regex, line) is not None:
                    # Search for associated KO.
                    ko_match = re.search(ko_regex, line).group(0)
                    if gene_id not in gene_to_kos:
                        gene_to_kos[gene_id] = [ko_match]
                    else:
                        gene_to_kos[gene_id].append(ko_match)

    return gene_to_kos


def create_kofamkoala_network(kofam_koala_folder, output_folder, reference_folder=False):
    """From a folder containing kofam koala result files, reconstruct draft metabolic networks.

    Args:
        kofam_koala_folder (str): path to kofam koala annotation files folder
        output_folder (str): path to the output folder
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|kofamkoala| Begin KEGG metabolism mapping for folder {0}.'.format(kofam_koala_folder))

    if reference_folder is not False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|kofamkoala| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|kofamkoala| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Download Uniprot metadata and create a json file containing them.
    options = {'kofam_koala_folder': kofam_koala_folder, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__
    options['tool_dependencies']['python_package']['networkx'] = networkx_version

    kegg2bipartitegraph_kofamkoala_metadata = {}
    kegg2bipartitegraph_kofamkoala_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_kofamkoala_metadata['reference_path'] = kegg_model_path

    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_json_hierarchy_path = os.path.join(kegg_model_path, 'kegg_hierarchy.json')
    kegg_json_model_path = os.path.join(kegg_model_path, 'kegg_metadata.json')

    # Read the reference KEGG sbml file.
    # Use it to create the organism sbml file.
    reader = libsbml.SBMLReader()
    reference_kegg_document = reader.readSBML(kegg_sbml_model_path)
    reference_kegg_model = reference_kegg_document.getModel()
    reference_model_fbc = reference_kegg_model.getPlugin('fbc')

    remove_gene_products = []
    for gene_product in reference_model_fbc.getListOfGeneProducts():
        remove_gene_products.append(gene_product.id)
    for gene_product in remove_gene_products:
        reference_model_fbc.removeGeneProduct(gene_product)

    reference_reactions = {reaction.id: reaction for reaction in reference_kegg_model.getListOfReactions()}
    reference_species = reference_kegg_model.getListOfSpecies()
    reference_groups = reference_kegg_model.getPlugin("groups")

    with open(kegg_json_model_path, 'r') as input_metadata_json:
        json_data = json.load(input_metadata_json)

    kegg2bipartitegraph_kofamkoala_metadata['json_reference_metadata'] = json_data

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

    # Create GRAPHML output folder.
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
    # Retrieve EC, KO and add reactions.
    for kofamkoala_result_file in os.listdir(kofam_koala_folder):
        taxon_reactions = {}

        base_name = os.path.splitext(kofamkoala_result_file)[0]
        kofamkoala_annotation_path = os.path.join(kofam_koala_folder, kofamkoala_result_file)
        gene_to_kos = read_kofam_koala_txt(kofamkoala_annotation_path)

        # Use KO found to be associated to reference protein to retrieve KEGG reaction.
        ko_added_reactions = []
        all_kos = []
        for gene_id in gene_to_kos:
            ko_ids = gene_to_kos[gene_id]
            all_kos.extend(ko_ids)

            for ko_id in ko_ids:
                if ko_id in ko_to_reactions:
                    reaction_ids = ko_to_reactions[ko_id]
                    for reaction_id in reaction_ids:
                        ko_added_reactions.append(reaction_id)
                        if reaction_id not in taxon_reactions:
                            taxon_reactions[reaction_id] = [gene_id]
                        else:
                            taxon_reactions[reaction_id].append(gene_id)
        logger.info('|kegg2bipartitegraph|kofamkoala| Added {0} reactions from {1} KO for organism {2} from Kofam Koala output.'.format(len(set(ko_added_reactions)), len(set(all_kos)) , base_name))

        total_added_reactions = list(taxon_reactions.keys())

        # Create pathway file contening pathway with reactions in the taxon.
        pathways_output_file_path = os.path.join(pathways_output_folder_path, base_name+'.tsv')
        organism_pathways = write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions)

        # Create module file contening module with reactions in the taxon.
        modules_output_file_path = os.path.join(modules_output_folder_path, base_name+'.tsv')
        organism_modules = write_module_file(kegg_modules, modules_output_file_path, total_added_reactions)
        class_data.append([base_name, *[len(class_representation[class_module].intersection(set(organism_modules))) for class_module in class_modules]])

        kegg_document, kegg_model = create_sbml_from_kegg_reactions(base_name, reference_reactions, reference_species, reference_groups, taxon_reactions, organism_pathways, organism_modules)

        stat_metabolic_networks[base_name] = (len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies()))

        # Create file if there is at least 1 reaction.
        if len(kegg_model.getListOfReactions()) > 0:
            # Create SBML file.
            sbml_output_file_path = os.path.join(sbml_output_folder_path, base_name+'.sbml')
            graphml_output_file_path = os.path.join(graphml_output_folder_path, base_name+'.graphml')
            libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
            sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
            logger.info('|kegg2bipartitegraph|kofamkoala| Network of {0} contains {1} reactions and {2} metabolites.'.format(base_name, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
        else:
            logger.info('|kegg2bipartitegraph|kofamkoala| No reactions in model for {0}, no SBML file will be created.'.format(base_name))

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

    kegg2bipartitegraph_kofamkoala_metadata['kegg2bipartitegraph_kofamkoala_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_kofamkoala_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_kofamkoala_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|kofamkoala| Draft networks creation complete.')
