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
import sys
import urllib.parse
import urllib.request
import libsbml

from networkx import __version__ as networkx_version
from bioservices import version as bioservices_version
from bioservices import KEGG

from kegg2bipartitegraph.utils import is_valid_dir, write_pathway_file, write_module_file
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def get_enzyme_org(organism):
    """ From a KEGG organism ID, retrieve the organism enzymes.

    Args:
        organism (str): KEGG code for an organism.

    Returns:
        enzymes (dict): Gene ID as key and list of Enzyme commission number as value.
    """
    response_text = KEGG_BIOSERVICES.link('enzyme', organism)
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    enzymes = {}
    for line in csvreader:
        gene_id = line[0].replace(organism+':', '')
        enzyme_id = line[1].replace('ec:', '')
        if gene_id not in enzymes:
            enzymes[gene_id] = [enzyme_id]
        else:
            enzymes[gene_id].append(enzyme_id)

    return enzymes


def get_ko_org(organism):
    """ From a KEGG organism ID, retrieve the organism KOs.

    Args:
        organism (str): KEGG code for an organism.

    Returns:
        enzymes (dict): Gene ID as key and list of KEGG Orthologs as value.
    """
    response_text = KEGG_BIOSERVICES.link('ko', organism)
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    kos = {}
    for line in csvreader:
        gene_id = line[0].replace(organism+':', '')
        ko_id = line[1].replace('ec:', '')
        if gene_id not in kos:
            kos[gene_id] = [ko_id]
        else:
            kos[gene_id].append(ko_id)

    return kos


def get_module_org(organism):
    """ From a KEGG organism ID, retrieve the organism modules.

    Args:
        organism (str): KEGG code for an organism.

    Returns:
        modules (dict): Gene ID as key and list of Module as value.
    """
    response_text = KEGG_BIOSERVICES.link('module', organism)
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    modules = {}
    for line in csvreader:
        gene_id = line[0].replace(organism+':', '')
        module_id = line[1].replace('md:'+organism+'_', '')
        if module_id not in modules:
            modules[module_id] = [gene_id]
        else:
            modules[module_id].append(gene_id)

    return modules


def get_pathway_org(organism):
    """ From a KEGG organism ID, retrieve the organism pathways.

    Args:
        organism (str): KEGG code for an organism.

    Returns:
        pathways (dict): Gene ID as key and list of pathway as value.
    """
    response_text = KEGG_BIOSERVICES.link('pathway', organism)
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    pathways = {}
    for line in csvreader:
        gene_id = line[0].replace(organism+':', '')
        pathway_id = line[1].replace('path:'+organism+'', 'map')
        if pathway_id not in pathways:
            pathways[pathway_id] = [gene_id]
        else:
            pathways[pathway_id].append(gene_id)

    return pathways


def create_organism_network(organism, output_folder, reference_folder=False):
    """From a KEGG organism ID create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number foundin organism and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        organism (str): KEGG code for an organism
        output_folder (str): path to the output folder
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|organism| Begin KEGG metabolism mapping for organism {0}.'.format(organism))

    if reference_folder != False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|organism| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|organism| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Create KEGG instance of bioservices.KEEG in this function to avoid trying to connect to KEGG with offline mode.
    global KEGG_BIOSERVICES
    KEGG_BIOSERVICES = KEGG()

    if KEGG_BIOSERVICES.isOrganism(organism) is not True:
        logger.info('|kegg2bipartitegraph|organism| Incorrect KEGG organism IDs {0}, please check at: https://www.genome.jp/kegg/catalog/org_list.html.'.format(organism))
        logger.info('|kegg2bipartitegraph|organism| Example: hsa for Homo sapiens (human).')
        sys.exit()
    # Download Uniprot metadata and create a json file containing them.
    options = {'organism': organism, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__
    options['tool_dependencies']['python_package']['networkx'] = networkx_version

    kegg2bipartitegraph_organism_metadata = {}
    kegg2bipartitegraph_organism_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_organism_metadata['reference_path'] = kegg_model_path

    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
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

    kegg2bipartitegraph_organism_metadata['json_reference_metadata'] = json_data

    # Retrieve organism pathway.
    org_kegg_pathways = get_pathway_org(organism)

    kegg_pathways = {}
    with open(kegg_pathways_path, 'r') as open_kegg_pathways_path:
        csvreader = csv.reader(open_kegg_pathways_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            pathway_id = line[0]
            pathway_name = line[1]
            pathway_reactions = line[2].split(',')
            if pathway_id in org_kegg_pathways:
                kegg_pathways[pathway_id] = (pathway_name, pathway_reactions)

    # Retrieve organism module.
    org_kegg_modules = get_module_org(organism)

    kegg_modules = {}
    with open(kegg_modules_path, 'r') as open_kegg_modules_path:
        csvreader = csv.reader(open_kegg_modules_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            module_id = line[0]
            module_name = line[1]
            module_reactions = line[3].split(',')
            if module_id in org_kegg_modules:
                kegg_modules[module_id] = (module_name, module_reactions)

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

    stat_metabolic_networks = {}
    # Retrieve proteins and annotations from esmecata annotation folder.
    enzymes = get_enzyme_org(organism)

    # Use EC found to be associated to reference protein to retrieve KEEG reaction.
    taxon_reactions = {}
    ec_added_reactions = []
    all_ecs = []
    for gene_id in enzymes:
        ec_ids = enzymes[gene_id]
        all_ecs.extend(ec_ids)

        for ec_id in ec_ids:
            if ec_id in ec_to_reactions:
                reaction_ids = ec_to_reactions[ec_id]
                for reaction_id in reaction_ids:
                    ec_added_reactions.append(reaction_id)
                    if reaction_id not in taxon_reactions:
                        taxon_reactions[reaction_id] = [gene_id]
                    else:
                        taxon_reactions[reaction_id].append(gene_id)
    logger.info('|kegg2bipartitegraph|organism| Added {0} reactions from {1} EC for organism {2}.'.format(len(set(ec_added_reactions)), len(set(all_ecs)), organism))

    total_added_reactions = list(taxon_reactions.keys())

    # Create pathway file contening pathway with reactions in the taxon.
    pathways_output_file_path = os.path.join(pathways_output_folder_path, organism+'.tsv')
    organism_pathways = write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions)

    # Create module file contening module with reactions in the taxon.
    modules_output_file_path = os.path.join(modules_output_folder_path, organism+'.tsv')
    organism_modules = write_module_file(kegg_modules, modules_output_file_path, total_added_reactions)
    
    kegg_document, kegg_model = create_sbml_from_kegg_reactions(organism, reference_reactions, reference_species, reference_groups, taxon_reactions, organism_pathways, organism_modules)

    stat_metabolic_networks[organism] = (len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies()))

    # Create file if there is at least 1 reaction.
    if len(kegg_model.getListOfReactions()) > 0:
        # Create SBML file.
        sbml_output_file_path = os.path.join(sbml_output_folder_path, organism+'.sbml')
        graphml_output_file_path = os.path.join(graphml_output_folder_path, organism+'.graphml')
        libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
        sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
        logger.info('|kegg2bipartitegraph|organism| Network of {0} contains {1} reactions and {2} metabolites.'.format(organism, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
    else:
        logger.info('|kegg2bipartitegraph|organism| No reactions in model for {0}, no SBML file will be created.'.format(organism))

    kegg_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    with open(kegg_stat_file, 'w') as stat_file_open:
        csvwriter = csv.writer(stat_file_open, delimiter='\t')
        csvwriter.writerow(['observation_name', 'Number_reactions', 'Number_metabolites'])
        for observation_name in stat_metabolic_networks:
            csvwriter.writerow([observation_name, stat_metabolic_networks[observation_name][0], stat_metabolic_networks[observation_name][1]])

    endtime = time.time()

    duration = endtime - starttime

    kegg2bipartitegraph_organism_metadata['kegg2bipartitegraph_organism_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_organism_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_organism_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|organism| Draft networks creation complete.')
