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
import logging
import os
import json
import time
import sys
import urllib.parse
import urllib.request
import libsbml

from bioservices import version as bioservices_version
from bioservices import KEGG

from kegg2bipartitegraph.utils import is_valid_dir
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)

# Create KEGG instance of bioservices.KEEG.
KEGG_BIOSERVICES = KEGG()

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def get_enzyme_org(organism):
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

def old_create_organism_sbml_cobrapy(organism, kegg_sbml_model_path, taxon_reactions, sbml_output_folder_path):
    """ Unused function only kept at this moment for archiving.
    """
    from cobra.io import read_sbml_model, write_sbml_model
    from cobra import Model
    kegg_model = read_sbml_model(kegg_sbml_model_path)

    species_model = Model(organism)

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
        sbml_output_file_path = os.path.join(sbml_output_folder_path, organism+'.sbml')
        write_sbml_model(species_model, sbml_output_file_path)
    else:
        logger.info('|kegg2bipartitegraph|organism| No reactions in model for {0}, no SBML file will be created.'.format(organism))


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


def create_organism_network(organism, output_folder):
    """From the output folder of 'esmecata annotation' create KEGG SBML files using bioservices.KEGG.
    To retrieve KEGG reactions, a mapping is performed between EC number and KEGG reactions.
    And if the option mapping_ko is set to True, it will also map KO ID to KEGG reaction

    Args:
        organism (str): KEGG code for an organism
        output_folder (str): path to the output folder
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|organism| Begin KEGG metabolism mapping for organism {0}.'.format(organism))

    if KEGG_BIOSERVICES.isOrganism(organism) is not True:
        logger.info('|kegg2bipartitegraph|organism| Incorrect KEGG organism IDs {0}, please check at: https://www.genome.jp/kegg/catalog/org_list.html.'.format(organism))
        logger.info('|kegg2bipartitegraph|organism| Example: hsa for Homo sapiens (human).')
        sys.exit()
    # Download Uniprot metadata and create a json file containing them.
    options = {'organism': organism, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['esmecata'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__

    kegg2bipartitegraph_organism_metadata = {}
    kegg2bipartitegraph_organism_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    data_kegg_model_path = DATA_ROOT
    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(data_kegg_model_path, 'kegg_model')
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_organism_metadata['reference_path'] = data_kegg_model_path

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
    modules_output_file_path = os.path.join(modules_output_folder_path, organism+'.tsv')
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
        sbml_output_file_path = os.path.join(sbml_output_folder_path, organism+'.sbml')
        graphml_output_file_path = os.path.join(graphml_output_folder_path, organism+'.graphml')
        libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
        sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
        logger.info('|kegg2bipartitegraph|organism| Network of {0} contains {1} reactions and {2} metabolites.'.format(organism, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
    else:
        logger.info('|kegg2bipartitegraph|organism| No reactions in model for {0}, no SBML file will be created.'.format(organism))


    clust_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    compute_stat_kegg(sbml_output_folder_path, clust_stat_file)

    endtime = time.time()

    duration = endtime - starttime

    kegg2bipartitegraph_organism_metadata['kegg2bipartitegraph_organism_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_organism_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_organism_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|organism| Draft networks creation complete.')
