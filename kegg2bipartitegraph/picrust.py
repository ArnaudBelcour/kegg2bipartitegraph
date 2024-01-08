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
import gzip
import libsbml
import time
import sys
import urllib.parse
import urllib.request

from networkx import __version__ as networkx_version
from bioservices import version as bioservices_version

from kegg2bipartitegraph.utils import is_valid_dir, write_pathway_file, write_module_file
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def create_picrust_network(input_folder, output_folder, reference_folder=False):
    """From the output folder of 'picrust_pipeline.py' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KO.

    Args:
        input_folder (str): path to the output folder of picrust2
        output_folder (str): path to the output folder
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|picrust| Begin KEGG metabolism mapping.')

    if reference_folder is not False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|picrust| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|picrust| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_folder': input_folder, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__
    options['tool_dependencies']['python_package']['networkx'] = networkx_version

    kegg2bipartitegraph_picrust_metadata = {}
    kegg2bipartitegraph_picrust_metadata['tool_options'] = options

    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_picrust_metadata['reference_path'] = kegg_model_path

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_json_model_path = os.path.join(kegg_model_path, 'kegg_metadata.json')
    kegg_json_hierarchy_path = os.path.join(kegg_model_path, 'kegg_hierarchy.json')

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

    kegg2bipartitegraph_picrust_metadata['json_reference_metadata'] = json_data

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

    # Retrieve EC and KO from picrust folder.
    compressed_ec_file = os.path.join(input_folder, 'EC_predicted.tsv.gz')
    ec_orgs = {}
    with gzip.open(compressed_ec_file,'rt') as open_ec_file:        
        csvreader = csv.reader(open_ec_file, delimiter='\t')
        header = next(csvreader)
        for line in csvreader:
            org_id = line[0]
            ec_numbers = [header[index+1].replace('EC:', '') for index, value in enumerate(line[1:]) if int(value) > 0]
            ec_orgs[org_id] = ec_numbers

    compressed_ko_file = os.path.join(input_folder, 'KO_predicted.tsv.gz')
    ko_orgs = {}
    with gzip.open(compressed_ko_file,'rt') as open_ec_file:        
        csvreader = csv.reader(open_ec_file, delimiter='\t')
        header = next(csvreader)
        for line in csvreader:
            org_id = line[0]
            ko_numbers = [header[index+1] for index, value in enumerate(line[1:]) if int(value) > 0]
            ko_orgs[org_id] = ko_numbers
    all_orgs = set(list(ec_orgs.keys()) + list(ko_orgs.keys())) 

    clust_pathways_output_folder_path = pathways_output_folder_path
    is_valid_dir(clust_pathways_output_folder_path)

    clust_sbml_output_folder_path = sbml_output_folder_path
    is_valid_dir(clust_sbml_output_folder_path)

    # Extract class modules.
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
    for org in all_orgs:
        taxon_reactions = {}

        # Map KO to KEGG reactions.
        ko_added_reactions = []
        for ko_id in ko_orgs[org]:
            if ko_id in ko_to_reactions:
                reaction_ids = ko_to_reactions[ko_id]
                for reaction_id in reaction_ids:
                    ko_added_reactions.append(reaction_id)
                    if reaction_id not in taxon_reactions:
                        taxon_reactions[reaction_id] = [ko_id]
                    else:
                        taxon_reactions[reaction_id].append(ko_id)
        logger.info('|kegg2bipartitegraph|picrust| Added {0} reactions from {1} KO for taxon {2}.'.format(len(set(ko_added_reactions)), len(ko_orgs[org]), org))

        # Use EC to retrieve KEEG reaction.
        ec_added_reactions = []
        for ec_id in ec_orgs[org]:
            if ec_id in ec_to_reactions:
                reaction_ids = ec_to_reactions[ec_id]
                for reaction_id in reaction_ids:
                    ec_added_reactions.append(reaction_id)
                    if reaction_id not in taxon_reactions:
                        taxon_reactions[reaction_id] = []
        logger.info('|kegg2bipartitegraph|picrust| Added {0} reactions from {1} EC for taxon {2}.'.format(len(set(ec_added_reactions)), len(ec_orgs[org]), org))

        total_added_reactions = list(taxon_reactions.keys())
        # Create pathway file contening pathway with reactions in the taxon.
        pathways_output_file_path = os.path.join(clust_pathways_output_folder_path, org+'.tsv')
        organism_pathways = write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions)

        # Create module file contening module with reactions in the taxon.
        modules_output_file_path = os.path.join(modules_output_folder_path, org+'.tsv')
        organism_modules = write_module_file(kegg_modules, modules_output_file_path, total_added_reactions)
        class_data.append([org, *[len(class_representation[class_module].intersection(set(organism_modules))) for class_module in class_modules]])

        kegg_document, kegg_model = create_sbml_from_kegg_reactions(org, reference_reactions, reference_species, reference_groups, taxon_reactions, organism_pathways, organism_modules)

        stat_metabolic_networks[org] = (len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies()))
        # Create file if there is at least 1 reaction.
        if len(kegg_model.getListOfReactions()) > 0:
            # Create SBML file.
            sbml_output_file_path = os.path.join(sbml_output_folder_path, org+'.sbml')
            libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
            graphml_output_file_path = os.path.join(graphml_output_folder_path, org+'.graphml')
            sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
            logger.info('|kegg2bipartitegraph|picrust| Network of {0} contains {1} reactions and {2} metabolites.'.format(org, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
        else:
            logger.info('|kegg2bipartitegraph|picrust| No reactions in model for {0}, no SBML file will be created.'.format(org))

    stat_number_kegg_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    with open(stat_number_kegg_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_reactions', 'Number_metabolites'])
        for observation_name in stat_metabolic_networks:
            csvwriter.writerow([observation_name, stat_metabolic_networks[observation_name][0], stat_metabolic_networks[observation_name][1]])

    module_class_file = os.path.join(output_folder, 'module_class.tsv')
    with open(module_class_file, 'w') as open_module_class_file:
        csvwriter = csv.writer(open_module_class_file, delimiter='\t')
        csvwriter.writerow(['observation_name', *class_modules])
        for class_dat in class_data:
            csvwriter.writerow(class_dat)

    endtime = time.time()

    duration = endtime - starttime
    kegg2bipartitegraph_picrust_metadata['kegg2bipartitegraph_picrust_duration'] = duration
    kegg2bipartitgraph_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_picrust_kegg.json')
    with open(kegg2bipartitgraph_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_picrust_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|picrust| Draft networks creation complete in {0}s.'.format(duration))
