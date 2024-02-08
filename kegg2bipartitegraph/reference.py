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
import json
import libsbml
import os
import re
import time
import sys
import urllib.parse
import urllib.request
import zipfile

from bioservices import version as bioservices_version
from bioservices import KEGG
from networkx import __version__ as networkx_version

from kegg2bipartitegraph.utils import is_valid_dir, urllib_query
from kegg2bipartitegraph.sbml import libsbml_check, initiate_sbml_model
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.graph import sbml_to_graphml

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])

# Ubiquitous metabolites defined in 10.1038/s41598-021-91486-8.
# 14 metabolites, consisting of: H2O, ATP, ADP, NAD+, NADH, NADP+, NADPH, CO2, ammonia, sulfate, thioredoxin, phosphate, pyrophosphate (PPi), and H+.
# They are removed by default from the graphs.
UBIQUITOUS_METABOLITES = ["C00001", "C00002", "C00003", "C00004", "C00005", "C00006",
                          "C00008", "C00009", "C00011", "C00013", "C00014", "C00059",
                          "C00080", "C00342"]


def get_reactions(reaction_folder):
    """Using bioservices.KEGG retrieve all keg files associated to KEGG reactions

    Args:
        reaction_folder (str): output folder which will contains reaction file
    """
    is_valid_dir(reaction_folder)

    response_text = KEGG_BIOSERVICES.list('reaction')

    reaction_ids = []
    for line in response_text.splitlines():
        reaction_id = line.split('\t')[0].strip()
        reaction_ids.append(reaction_id)

    already_downloaded_reaction = [kegg_file.replace('.keg', '') for kegg_file in os.listdir(reaction_folder)]
    reaction_ids = sorted(list(set(reaction_ids) - set(already_downloaded_reaction)))
    logger.info('|kegg2bipartitegraph|reference| Retrieving {0} reactions from KEGG.'.format(len(reaction_ids)))
    # Download reaction by chunk of 9.
    for i in range(0, len(reaction_ids), 9):
        reaction_chunk = '+'.join(reaction_ids[i:i+9])
        response_text = KEGG_BIOSERVICES.get(reaction_chunk)
        if isinstance(response_text, int):
            time.sleep(600)
            response_text = KEGG_BIOSERVICES.get(reaction_chunk)
        for response_split in response_text.split('///')[:-1]:
            if response_split.startswith('\n'):
                response_split = response_split.replace('\n', '', 1)
            reaction_data = KEGG_BIOSERVICES.parse(response_split.strip('\n'))

            reaction_id = reaction_data['ENTRY'].split(' ')[0]
            reaction_file = os.path.join(reaction_folder, reaction_id+'.keg')
            if not os.path.exists(reaction_file):
                with open(reaction_file, 'w') as output_file:
                    output_file.write(response_split+'///')
                time.sleep(3)


def extract_reaction(reaction_id, equation_text):
    """Using the EQUATION in the reaction keg file, extract the compounds from the reaction formula.

    Args:
        reaction_id (str): Id of the reaction
        equation_text (str): string containing the reaction formula

    Returns:
        left_compounds (list): compounds at the left of the reaction formula
        right_compounds (list): compounds at the right of the reaction formula
        modify_stoichiometry (bool): if stoichiometry of reaction changed
    """
    if '<=>' not in equation_text:
        logger.critical('|kegg2bipartitegraph|reference| No <=> symbol in equation {0} of {1}.'.format(equation_text, reaction_id))
    if len(equation_text.split('<=>')) > 2:
        logger.critical('|kegg2bipartitegraph|reference| More than one <=> symbol in equation {0} of {1}.'.format(equation_text, reaction_id))

    # Find group pattern
    equation_pattern = r'(?P<stoichiometry>\d|\w|\([\w\d\+]*\))*\ *(?P<compound>[CG]\d{5})|(?P<symbol><*=>*)'
    left_compounds = []
    right_compounds = []
    equation_symbol_passed = False
    modify_stoichiometry = False
    for m in re.finditer(equation_pattern, equation_text):
        stoichiometry = m.groupdict()['stoichiometry']
        # If no stoichiometry is indicated, it is a stoichiometry of 1.
        if stoichiometry is None:
            stoichiometry = 1
        else:
            try:
                stoichiometry = int(stoichiometry)
            except:
                logger.critical('|kegg2bipartitegraph|reference| Stoichiometry is not an int ({0}) for {1}, replace by 1 as it is not relevant for topological analysis.'.format(stoichiometry, reaction_id))
                stoichiometry = 1
                modify_stoichiometry = True
        compound = m.groupdict()['compound']
        symbol = m.groupdict()['symbol']
        if symbol is not None:
            equation_symbol_passed = True
        if equation_symbol_passed is False and compound is not None:
            left_compounds.append((compound, stoichiometry))
        elif equation_symbol_passed is True and compound is not None:
            right_compounds.append((compound, stoichiometry))

    return left_compounds, right_compounds, modify_stoichiometry


def get_compound_names(compound_file):
    """Using bioservices.KEGG find the name associated to compound and glycan.

    Args:
        compound_file (str): output file which will contains compound ID and name
    """
    response_text = KEGG_BIOSERVICES.list('compound')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    compounds = {}
    for line in csvreader:
        cpd_id = line[0].replace('cpd:', '')
        cpd_name = line[1].split('; ')[0]
        compounds[cpd_id] = cpd_name

    response_text = KEGG_BIOSERVICES.list('glycan')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []

    for line in csvreader:
        cpd_id = line[0].replace('gl:', '')
        cpd_name = line[1].split('; ')[0]
        compounds[cpd_id] = cpd_name

    with open(compound_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['compound_id', 'compound_name'])
        for cpd_id in compounds:
            csvwriter.writerow([cpd_id, compounds[cpd_id]])


def get_matching_elements(database_1, database_2, replace_id=None):
    """Query KEGG suing link to match 2 databases (for example module and reaction) and returns a dict with matching element.
    For example get_matching_elements('module', 'ko'): element_1_element_2['module_1'] = ['ko_1', 'ko_2']

    Args:
        database_1 (str): name of first database to match
        database_2 (str): name of second database to match
        replace_id (list): str to replace in element_id_1 (first element of list replace by second)
    Returns:
        element_1_element_2 (dict): dictionary of matching with element from database 1 as key and list of element of database 2 as key.
    """
    response_text = KEGG_BIOSERVICES.link(database_1, database_2)
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    element_1_element_2 = {}
    for line in csvreader:
        element_2_id = line[0].split(':')[1]
        if replace_id is None:
            element_1_id = line[1].split(':')[1]
        else:
            element_1_id = line[1].split(':')[1].replace(replace_id[0], replace_id[1])
        if element_1_id not in element_1_element_2:
            element_1_element_2[element_1_id] = [element_2_id]
        else:
            element_1_element_2[element_1_id].append(element_2_id)

    return element_1_element_2


def get_modules(module_file):
    """Using bioservices.KEGG to find the module associated with reactions.

    Args:
        module_file (str): output file which will contains module ID, name, formula and reactions
    Returns:
        module_data (dict): dictionary with module ID as key and module name, formula, reactions, kos and compounds as key.
    """
    module_reactions = get_matching_elements('module', 'reaction')
    module_kos = get_matching_elements('module', 'ko')
    module_compounds = get_matching_elements('module', 'compound')

    response_text = KEGG_BIOSERVICES.list('module')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []

    modules = {}
    for line in csvreader:
        module_id = line[0]
        module_name = line[1]
        if '=>' in line[1]:
            module_data = line[1].split(', ')
            module_formula = module_data[-1]
        else:
            module_formula = ''
        modules[module_id] = (module_name, module_formula)

    all_modules = set(list(modules.keys()) + list(module_reactions.keys()))
    module_data= {}
    with open(module_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['module_id', 'module_name', 'module_formula', 'module_reactions', 'module_kos', 'module_compounds'])
        for module_id in all_modules:
            if module_id in module_reactions:
                module_reaction = ','.join(module_reactions[module_id])
            else:
                module_reaction = ''
            if module_id in module_kos:
                module_ko = ','.join(module_kos[module_id])
            else:
                module_ko = ''
            if module_id in module_compounds:
                module_compound = ','.join(module_compounds[module_id])
            else:
                module_compound = ''
            if module_id in modules:
                module_name = modules[module_id][0]
                module_formula = modules[module_id][1]
            else:
                module_name = ''
                module_formula = ''
            csvwriter.writerow([module_id, module_name, module_formula, module_reaction, module_ko, module_compound])
            module_data[module_id] = [module_name, module_formula, module_reaction, module_ko, module_compound]

    return module_data


def get_kegg_hierarchy(hierarchy_file):
    """Using bioservices.KEGG to create hierarchy of metabolite IDs.

    Args:
        hierarchy_file (str): output file which will contains hierarchy of KEGG IDs
    """
    # Create hierarchy for module ID.
    response_text = KEGG_BIOSERVICES.get('br:ko00002')

    module_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            module_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            module_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            module_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            module_id = d_category.split('  ')[0]
            module_hierarchy[a_category][b_category][c_category].append(module_id)
    time.sleep(3)

    # Create hierarchy for pathway.
    response_text = KEGG_BIOSERVICES.get('br:br08901')
    pathway_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            pathway_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            pathway_hierarchy[a_category][b_category] = []
        if line.startswith('C    '):
            pathway_id = 'map' + line.replace('C    ', '').split('  ')[0]
            pathway_hierarchy[a_category][b_category].append(pathway_id)
    time.sleep(3)

    # Create hierarchy for KO hierarchy.
    response_text = KEGG_BIOSERVICES.get('br:ko00001')

    ko_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            ko_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            ko_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            ko_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            ko_id = d_category.split('  ')[0]
            ko_hierarchy[a_category][b_category][c_category].append(ko_id)
    time.sleep(3)

    # Create hierarchy for KO transporter.
    response_text = KEGG_BIOSERVICES.get('br:ko02000')

    transporter_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            transporter_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            transporter_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            transporter_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            ko_id = d_category.split('  ')[0]
            transporter_hierarchy[a_category][b_category][c_category].append(module_id)
    time.sleep(3)

    # Create hierarchy for metabolite.
    response_text = KEGG_BIOSERVICES.get('br:br08001')

    metabolite_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            metabolite_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            metabolite_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            metabolite_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            metabolite_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # Create hierarchy for lipids.
    response_text = KEGG_BIOSERVICES.get('br:br08002')

    lipid_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            lipid_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            lipid_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            lipid_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            lipid_id = d_category.split('  ')[0]
            lipid_hierarchy[a_category][b_category][c_category].append(lipid_id)
    time.sleep(3)

    global_hierarchy = {}
    global_hierarchy['module'] = module_hierarchy
    global_hierarchy['pathway'] = pathway_hierarchy
    global_hierarchy['ko'] = ko_hierarchy
    global_hierarchy['ko_transporter'] = transporter_hierarchy
    global_hierarchy['metabolite'] = metabolite_hierarchy
    global_hierarchy['metabolite_lipid'] = lipid_hierarchy

    with open(hierarchy_file, 'w') as open_json_file:
        json.dump(global_hierarchy, open_json_file, indent=4)


def get_pathways(pathway_file):
    """Using bioservices.KEGG to find the pathway associated with reactions.

    Args:
        pathway_file (str): output file which will contains pathway ID, name, formula and reactions.
    Returns:
        pathway_data (dict): dictionary with pathway ID as key and pathway name, reactions, kos and compounds as key.
    """
    # Get pathway matching with reactions.
    pathway_reactions = get_matching_elements('pathway', 'reaction', ['rn', 'map'])
    # Get pathway matching with KOs.
    pathway_kos = get_matching_elements('pathway', 'ko', ['rn', 'map'])
    # Get pathway matching with compounds.
    pathway_compounds = get_matching_elements('pathway', 'compound', ['rn', 'map'])

    # Get pathway names.
    response_text = KEGG_BIOSERVICES.list('pathway')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    pathways = {}
    for line in csvreader:
        pathway_id = line[0]
        pathway_name = line[1]
        pathways[pathway_id] = pathway_name

    # Write pathway file.
    all_pathways = set(list(pathways.keys()) + list(pathway_reactions.keys()))
    pathway_data= {}
    with open(pathway_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_reactions', 'pathway_kos', 'pathway_compounds'])
        for pathway_id in all_pathways:
            if pathway_id in pathway_reactions:
                pathway_reaction = ','.join(pathway_reactions[pathway_id])
            else:
                pathway_reaction = ''
            if pathway_id in pathway_kos:
                pathway_ko = ','.join(pathway_kos[pathway_id])
            else:
                pathway_ko = ''
            if pathway_id in pathway_compounds:
                pathway_cpd = ','.join(pathway_compounds[pathway_id])
            else:
                pathway_cpd = ''
            if pathway_id in pathways:
                pathway_name = pathways[pathway_id]
            else:
                pathway_name = ''
            csvwriter.writerow([pathway_id, pathway_name, pathway_reaction, pathway_ko, pathway_cpd])
            pathway_data[pathway_id] = [pathway_name, pathway_reaction, pathway_ko, pathway_cpd]

    return pathway_data


def create_sbml_model_from_kegg_file_libsbml(reaction_folder, compound_file, output_sbml, output_tsv,
                                             kegg_removed_changed_reaction_path, pathway_data, module_data,
                                             remove_ubiquitous=True, remove_glycan_reactions=True):
    """Using the reaction keg files (from retrieve_reactions), the compound file (from get_compound_names),
    create a SBML file (containing all reactions of KEGG) and a tsv file (used to map EC and/or KO to reaction ID)
    Use python-libsbml instead of cobrapy as cobrapy does not handle metabolite being both in reactant and product of reaction:
    https://github.com/opencobra/cobrapy/issues/906#issuecomment-527686017

    Args:
        reaction_folder (str): path to the folder containing keg reaction files
        compound_file (str): path to the tsv file containing compound ID and name
        output_sbml (str): path to the sbml output file
        output_tsv (str): path to an output tsv file mapping reaction ID, with KO and EC
        kegg_removed_changed_reaction_path (str): path to an output tsv showing removed or modified reaction
        pathway_data (dict): dictionary with pathway ID as key, list of pathway name and reaction as values
        module_data (dict): dictionary with module ID as key, list of module name, formula and reaction as values
        remove_ubiquitous (bool): remove ubiquitous metabolites
        remove_glycan_reactions (bool): remove glycan associated reactions
    """
    # Initiate libsbml model.
    document, model, model_fbc, model_groups = initiate_sbml_model('KEGG')

    # Read compounds file.
    compounds = {}
    with open(compound_file, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compounds[line[0]] = line[1]

    reaction_ecs = {}
    reactions = {}
    genes = []
    pathways = {}
    metabolites = []
    remove_or_modify_reactions = {}

    # Parse reaction file to extract information.
    for reaction_file in os.listdir(reaction_folder):
        reaction_file_path = os.path.join(reaction_folder, reaction_file)
        with open(reaction_file_path) as open_reaction_file_path:
            reaction_data = KEGG_BIOSERVICES.parse(open_reaction_file_path.read())

        reaction_id = reaction_data['ENTRY'].split(' ')[0]
        if 'NAME' in reaction_data:
            reaction_name = reaction_data['NAME'][0]
        else:
            reaction_name = None
        left_compounds, right_compounds, modify_stoichiometry = extract_reaction(reaction_id, reaction_data['EQUATION'])
        if 'ORTHOLOGY' in reaction_data:
            kegg_orthologs = reaction_data['ORTHOLOGY'].keys()
        else:
            kegg_orthologs = []
        if 'ENZYME' in reaction_data:
            kegg_ecs = reaction_data['ENZYME']
        else:
            kegg_ecs = []
        if 'PATHWAY' in reaction_data:
            kegg_pathways = reaction_data['PATHWAY']
            for pathway in kegg_pathways:
                if pathway not in pathways:
                    pathways[pathway] = {}
                    pathways[pathway]['name'] = kegg_pathways[pathway]
                    pathways[pathway]['reaction'] = [reaction_id]
                else:
                    pathways[pathway]['reaction'].append(reaction_id)
        reaction_ecs[reaction_id] = [kegg_orthologs, kegg_ecs]
        if left_compounds is None:
            continue

        reactions[reaction_id] = {}
        glycan_reaction = None

        # Create SBML reaction.
        reaction = model.createReaction()
        libsbml_check(reaction, 'create reaction')
        libsbml_check(reaction.setId(reaction_id), 'set reaction id %s' %reaction_id)
        libsbml_check(reaction.setFast(False), 'set fast')
        reaction.setReversible(False)
        # Set FBC plugin for reaction.
        r_fbc: "libsbml.FbcReactionPlugin" = reaction.getPlugin("fbc")
        r_fbc.setLowerFluxBound('default_zero_bound')
        r_fbc.setUpperFluxBound('default_upper_bound')

        reactants = []
        # Create metabolites from left compounds, remove ubiquitous and mark glycan reactions.
        for stoichiometry_metabolite in left_compounds:
            metabolite_id = stoichiometry_metabolite[0]
            if metabolite_id.startswith('G'):
                glycan_reaction = True
                if remove_glycan_reactions is True:
                    continue
            if remove_ubiquitous is True:
                if metabolite_id in UBIQUITOUS_METABOLITES:
                    continue
            if metabolite_id not in metabolites:
                # Create species if not in SBML.
                s = model.createSpecies()
                libsbml_check(s, 'create species')
                libsbml_check(s.setId(metabolite_id), 'set species id %s' %metabolite_id)
                libsbml_check(s.setMetaId(metabolite_id), 'set species meta id %s' %metabolite_id)
                libsbml_check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
                libsbml_check(s.setHasOnlySubstanceUnits(False), 'set setHasOnlySubstanceUnits to False')
                libsbml_check(s.setConstant(False), 'set setConstant to False')
                libsbml_check(s.setInitialAmount(0.0), 'set initAmount')
                libsbml_check(s.setName(compounds[metabolite_id]), 'set species Name {0}'.format(compounds[metabolite_id]))
                libsbml_check(s.setCompartment('c'), 'set species compartment c')
                metabolites.append(metabolite_id)
            # Create species as reactant in SBML.
            species_ref = reaction.createReactant()
            libsbml_check(species_ref, 'create reactant')
            libsbml_check(species_ref.setSpecies(metabolite_id), 'assign reactant species %s' %metabolite_id)
            libsbml_check(species_ref.setStoichiometry(stoichiometry_metabolite[1]), 'set stoichiometry {0}'.format(stoichiometry_metabolite[1]))
            libsbml_check(species_ref.setConstant(False), 'set constant %s' %False)
            reactants.append(metabolite_id)

        products = []
        # Create metabolites from right compounds, remove ubiquitous and mark glycan reactions.
        for stoichiometry_metabolite in right_compounds:
            metabolite_id = stoichiometry_metabolite[0]
            if metabolite_id.startswith('G'):
                glycan_reaction = True
                if remove_glycan_reactions is True:
                    continue
            if remove_ubiquitous is True:
                if metabolite_id in UBIQUITOUS_METABOLITES:
                    continue
            if metabolite_id not in metabolites:
                # Create species if not in SBML.
                s = model.createSpecies()
                libsbml_check(s, 'create species')
                libsbml_check(s.setId(metabolite_id), 'set species id %s' %metabolite_id)
                libsbml_check(s.setMetaId(metabolite_id), 'set species meta id %s' %metabolite_id)
                libsbml_check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
                libsbml_check(s.setHasOnlySubstanceUnits(False), 'set setHasOnlySubstanceUnits to False')
                libsbml_check(s.setConstant(False), 'set setConstant to False')
                libsbml_check(s.setInitialAmount(0.0), 'set initAmount')
                libsbml_check(s.setName(compounds[metabolite_id]), 'set species Name {0}'.format(compounds[metabolite_id]))
                libsbml_check(s.setCompartment('c'), 'set species compartment c')
                metabolites.append(metabolite_id)
            # Create species as product in SBML.
            species_ref = reaction.createProduct()
            libsbml_check(species_ref, 'create product')
            libsbml_check(species_ref.setSpecies(metabolite_id), 'assign product species %s' %metabolite_id)
            libsbml_check(species_ref.setStoichiometry(stoichiometry_metabolite[1]), 'set stoichiometry {0}'.format(stoichiometry_metabolite[1]))
            libsbml_check(species_ref.setConstant(False), 'set constant %s' %False)
            products.append(metabolite_id)

        if reaction_name:
            reaction.name = reaction_name
        if kegg_orthologs != []:
            gpr_association = r_fbc.createGeneProductAssociation()

            for ko in kegg_orthologs:
                if ko not in genes:
                    # Create gene in gene product list of SBML.
                    gene_prod = model_fbc.createGeneProduct()
                    libsbml_check(gene_prod.setId(ko), 'add gene %s' %ko)
                    gene_prod.setName(ko)
                    gene_prod.setLabel(ko)
                    genes.append(ko)

            # Add gene to reaction with GPR.
            gpr = ' or '.join([ko for ko in kegg_orthologs])
            libsbml_check(gpr_association.setAssociation(gpr, True, True), "set gpr: " + gpr)
        else:
            gpr = ''

        remove_reaction = False
        reactions_metabolites = reactants + products
        # Remove reactions without reactants and products.
        if reactions_metabolites == []:
            logger.critical('|kegg2bipartitegraph|reference| No reactants and products for {0}, will be removed from model.'.format(reaction_id))
            remove_or_modify_reactions[reaction_id] = ['no_reactant_and_product', reaction_name, ','.join(reactants), ','.join(products), left_compounds, right_compounds, gpr]
            remove_reaction = True
        # Remvoe glycan reactions.
        if remove_glycan_reactions is True and glycan_reaction is True:
            logger.critical('|kegg2bipartitegraph|reference| Do not add glycan reaction {0}.'.format(reaction_id))
            remove_or_modify_reactions[reaction_id] = ['glycan', reaction_name, ','.join(reactants), ','.join(products), left_compounds, right_compounds, gpr]
            remove_reaction = True
        # Remove reactions if all left compounds are UBIQUITOUS_METABOLITES, because it unlocks metabolites freely.
        if all([compound[0] in UBIQUITOUS_METABOLITES for compound in left_compounds]):
            logger.critical('|kegg2bipartitegraph|reference| Remove reaction {0} as all its reactants are ubiquitous metabolites.'.format(reaction_id))
            remove_or_modify_reactions[reaction_id] = ['ubiquitous_reactants', reaction_name, ','.join(reactants), ','.join(products), left_compounds, right_compounds, gpr]
            remove_reaction = True
        # Remove reactions if all right compounds are UBIQUITOUS_METABOLITES, because it unlocks metabolites freely.
        if all([compound[0] in UBIQUITOUS_METABOLITES for compound in right_compounds]):
            logger.critical('|kegg2bipartitegraph|reference| Remove reaction {0} as all its products are ubiquitous metabolites.'.format(reaction_id))
            remove_or_modify_reactions[reaction_id] = ['ubiquitous_products', reaction_name, ','.join(reactants), ','.join(products), left_compounds, right_compounds, gpr]
            remove_reaction = True

        # Remove reaction if it contains glycan metabolite or it does not have reactants and products.
        if remove_reaction is True:
            model.removeReaction(reaction_id)

        if modify_stoichiometry:
            remove_or_modify_reactions[reaction_id] = ['change_stoichiometry', reaction_name, ','.join(reactants), ','.join(products), left_compounds, right_compounds, gpr]

    # Add pathways using groups:kind.
    # Help from: https://synonym.caltech.edu/software/libsbml/5.18.0/docs/formatted/python-api/groups_example1_8py-example.html
    for pathway_id in pathway_data:
        pathway_name = pathway_data[pathway_id][0]
        pathway_reactions = pathway_data[pathway_id][1].split(',')
        if pathway_reactions != ['']:
            group = model_groups.createGroup()
            group.setId(pathway_id)
            group.setName(pathway_name)
            group.setKind("partonomy")
            for reaction in pathway_reactions:
                member = group.createMember()
                member.setId(reaction)
                member.setIdRef(reaction)

    # Add module using groups:kind.
    for module_id in module_data:
        module_name = module_data[module_id][0]
        module_reactions = module_data[module_id][2].split(',')
        if module_reactions != ['']:
            group = model_groups.createGroup()
            group.setId(module_id)
            group.setName(module_name)
            group.setKind("partonomy")
            for reaction in module_reactions:
                member = group.createMember()
                member.setId(reaction)
                member.setIdRef(reaction)

    libsbml.writeSBMLToFile(document, output_sbml)
    logger.info('|kegg2bipartitegraph|reference| {0} reactions and {1} metabolites in reference model.'.format(len(model.getListOfReactions()), len(model.getListOfSpecies())))
    removed_reactions = set([rxn_id for rxn_id in remove_or_modify_reactions if remove_or_modify_reactions[rxn_id][0] != 'change_stoichiometry'])
    logger.info('|kegg2bipartitegraph|reference| Remove {0} reactions in reference model (cause: glycan, ubiquitous metabolites, no reactants and products).'.format(len(removed_reactions)))
    modified_reactions = set([rxn_id for rxn_id in remove_or_modify_reactions if remove_or_modify_reactions[rxn_id][0] == 'change_stoichiometry'])
    logger.info('|kegg2bipartitegraph|reference| Modify {0} reaction stoichiometry in reference model.'.format(len(modified_reactions)))

    with open(output_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['reaction_id', 'kegg_orthologs', 'kegg_ec'])
        for reaction_id in reaction_ecs:
            csvwriter.writerow([reaction_id, ','.join(reaction_ecs[reaction_id][0]), ','.join(reaction_ecs[reaction_id][1])])

    with open(kegg_removed_changed_reaction_path, 'w') as open_kegg_removed_changed_reaction_path:
        csvwriter = csv.writer(open_kegg_removed_changed_reaction_path, delimiter='\t')
        csvwriter.writerow(['reaction_id', 'reason', 'reaction_name', 'reactants', 'products', 'stoichiometry_left', 'stoichiometry_right', 'GPR'])
        for reaction_id in remove_or_modify_reactions:
            csvwriter.writerow([reaction_id, *remove_or_modify_reactions[reaction_id]])


def get_kegg_database_version():
    """Retrieve the KEGG version release

    Returns:
        kegg_version (str): string containing KEGG release version
    """
    try:
        KEGG_BIOSERVICES
    except NameError:
        KEGG_BIOSERVICES = KEGG()

    response_text = KEGG_BIOSERVICES.dbinfo()

    for line in response_text.splitlines():
        if 'Release ' in line:
            kegg_version = line.split('Release ')[1].strip()

    return kegg_version


def go_to_ec(ec_to_gos_file):
    """Download the ec2go and extracts the association between EC and GO into a file.

    Args:
        ec_to_gos_file (str): output file which will contains association between EC and GO terms
    """
    regex = r'GO:\d{7}'

    response = urllib_query('http://current.geneontology.org/ontology/external2go/ec2go')
    ec2gos = {}
    for line in response.readlines():
        read_line = line.decode('utf-8')
        if not read_line.startswith('!'):

            ec, gos = read_line.split(' >')
            ec = ec.replace('EC:', '')
            if ec in ec2gos:
                print(ec)
            go_id = re.search(regex, gos).group(0)
            if '-' not in ec:
                if ec not in ec2gos:
                    ec2gos[ec] = [go_id]
                else:
                    ec2gos[ec].append(go_id)

    output_tsv = os.path.join(ec_to_gos_file)
    with open(output_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['ec_id', 'go_id'])
        for ec_id in ec2gos:
            csvwriter.writerow([ec_id, ','.join(ec2gos[ec_id])])


def create_reference_base(output_folder=None):
    """ Create kegg2bipartiegraph database by querying KEGG API and processing it.

    Args:
        output_folder (str): path to the folder that will contain reference data, by default it is in the package repository.
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|reference| Begin KEGG metabolism reference model creation.')

    # Create KEGG instance of bioservices.KEEG in this function to avoid trying to connect to KEGG with offline mode.
    global KEGG_BIOSERVICES
    KEGG_BIOSERVICES = KEGG()

    # Create a json file containing metadata.
    options = {}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__
    options['tool_dependencies']['python_package']['networkx'] = networkx_version

    kegg2bipartitegraph_reference_metadata = {}
    kegg2bipartitegraph_reference_metadata['tool_options'] = options
    kegg2bipartitegraph_reference_metadata['kegg_release_number'] = get_kegg_database_version()

    # If no path is indicated, create reference data in package folder.
    if output_folder is None:
        output_folder = DATA_ROOT

    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(output_folder, 'kegg_model')
    is_valid_dir(kegg_model_path)

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    kegg_compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_graphml_model_path = os.path.join(kegg_model_path, 'kegg_model.graphml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_removed_changed_reaction_path = os.path.join(kegg_model_path, 'kegg_removed_changed_reaction.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_hierarchy_path = os.path.join(kegg_model_path, 'kegg_hierarchy.json')
    kegg_metadata_path = os.path.join(kegg_model_path, 'kegg_metadata.json')
    ec2gos_file = os.path.join(kegg_model_path, 'ec_to_gos.tsv')

    logger.info('|kegg2bipartitegraph|reference| Check missing files in {0}.'.format(DATA_ROOT))
    input_files = [kegg_compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path,
                   kegg_pathways_path, kegg_modules_path, kegg_reactions_folder_path,
                   kegg_removed_changed_reaction_path, kegg_hierarchy_path,
                   ec2gos_file]

    missing_files = []
    for input_file in input_files:
        if not os.path.exists(input_file):
            missing_files.append(input_file)
    
    if len(missing_files) > 0:
        logger.info('|kegg2bipartitegraph|reference| Missing: ' + ' '.join(missing_files))
        if not os.path.exists(kegg_reactions_folder_path):
            logger.info('|kegg2bipartitegraph|reference| Retrieve reactions from KEGG to create SBML model.')
            get_reactions(kegg_reactions_folder_path)
        else:
            if len(os.listdir(kegg_reactions_folder_path)) != len(KEGG_BIOSERVICES.reactionIds):
                logger.info('|kegg2bipartitegraph|reference| Retrieve reactions from KEGG to create SBML model.')
                get_reactions(kegg_reactions_folder_path) 

        if not os.path.exists(ec2gos_file):
            go_to_ec(ec2gos_file)
        if not os.path.exists(kegg_compound_file_path):
            logger.info('|kegg2bipartitegraph|reference| Retrieve compound IDs and names from KEGG to create SBML model.')
            get_compound_names(kegg_compound_file_path)

        if not os.path.exists(kegg_modules_path):
            logger.info('|kegg2bipartitegraph|reference| Create KEGG reference module file.')
            module_data = get_modules(kegg_modules_path)
        else:
            module_data = {}
            with open(kegg_modules_path, 'r') as open_kegg_modules_path:
                csvreader = csv.reader(open_kegg_modules_path, delimiter='\t')
                next(csvreader)
                for line in csvreader:
                    module_data[line[0]] = line[1:]

        if not os.path.exists(kegg_hierarchy_path):
            logger.info('|kegg2bipartitegraph|reference| Create KEGG reference module hierarchy file.')
            get_kegg_hierarchy(kegg_hierarchy_path)

        if not os.path.exists(kegg_pathways_path):
            logger.info('|kegg2bipartitegraph|reference| Create KEGG reference pathway file.')
            pathway_data = get_pathways(kegg_pathways_path)
        else:
            pathway_data = {}
            with open(kegg_pathways_path, 'r') as open_kegg_pathways_path:
                csvreader = csv.reader(open_kegg_pathways_path, delimiter='\t')
                next(csvreader)
                for line in csvreader:
                    pathway_data[line[0]] = line[1:]

        logger.info('|kegg2bipartitegraph|reference| Create KEGG reference SBML and mapping tsv file.')
        create_sbml_model_from_kegg_file_libsbml(kegg_reactions_folder_path, kegg_compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path, kegg_removed_changed_reaction_path,
                                                 pathway_data, module_data)
        sbml_to_graphml(kegg_sbml_model_path, kegg_graphml_model_path)

        # Create compress archive.
        model_zipfile = zipfile.ZipFile(KEGG_ARCHIVE, mode="w")
        for filepath in input_files:
            file_name = os.path.basename(filepath)
            model_zipfile.write(filepath, file_name)
        model_zipfile.close()
    else:
        logger.info('|kegg2bipartitegraph|reference| No missing files.')

    endtime = time.time()

    duration = endtime - starttime

    kegg2bipartitegraph_reference_metadata['kegg2bipartitegraph_reference_duration'] = duration
    with open(kegg_metadata_path, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_reference_metadata, ouput_file, indent=4)

    logger.info('|kegg2bipartitegraph|reference| Reference creation finished in {0}.'.format(duration))
