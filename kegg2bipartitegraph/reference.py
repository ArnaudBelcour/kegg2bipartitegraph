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
import re
import time
import sys
import urllib.parse
import urllib.request
import zipfile

from bioservices import version as bioservices_version
from bioservices import KEGG
from cobra import __version__ as cobra_version
from cobra import Model, Reaction, Metabolite
from cobra.io import write_sbml_model

from kegg2bipartitegraph.utils import is_valid_dir
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)
logging.getLogger("cobra.io.sbml").setLevel(logging.CRITICAL)

# Create KEGG instance of bioservices.KEEG.
KEGG_BIOSERVICES = KEGG()

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

    for reaction_id in reaction_ids:
        response_text = KEGG_BIOSERVICES.get(reaction_id)
        reaction_file = os.path.join(reaction_folder, reaction_id+'.keg')
        with open(reaction_file, 'w') as output_file:
            output_file.write(response_text)


def extract_reaction(reaction_id, equation_text):
    """Using the EQUATION in the reaction keg file, extract the compounds from the reaction formula.

    Args:
        reaction_id (str): Id of the reaction
        equation_text (str): string containing the reaction formula

    Returns:
        left_compounds (list): compounds at the left of the reaction formula
        right_compounds (list): compounds at the right of the reaction formula
    """
    equation_pattern = r'(?P<stoechiometry>\d|\w|\([\w\d\+]*\))*\ *(?P<compound>[CG]\d{5})|(?P<symbol><*=>*)'
    left_compounds = []
    right_compounds = []
    equation_symbol_passed = False
    for m in re.finditer(equation_pattern, equation_text):
        stoechiometry = m.groupdict()['stoechiometry']
        if stoechiometry is None:
            stoechiometry = 1
        else:
            try:
                stoechiometry = int(stoechiometry)
            except:
                logger.critical('Stochiometry is not an int ({0}) for {1}, replace by 1 as it is not relevant for topological analysis.'.format(stoechiometry, reaction_id))
                stoechiometry = 1
        compound = m.groupdict()['compound']
        symbol = m.groupdict()['symbol']
        if symbol is not None:
            equation_symbol_passed = True
        if equation_symbol_passed is False and compound is not None:
            left_compounds.append((compound, stoechiometry))
        elif equation_symbol_passed is True and compound is not None:
            right_compounds.append((compound, stoechiometry))

    return left_compounds, right_compounds


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


def get_modules(module_file):
    """Using bioservices.KEGG to find the module associated with reacitons.

    Args:
        module_file (str): output file which will contains module ID, name, formula and reactions
    """
    response_text = KEGG_BIOSERVICES.link('module', 'reaction')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []
    module_reactions = {}
    for line in csvreader:
        reaction_id = line[0].split(':')[1]
        module_id = line[1].split(':')[1]
        if module_id not in module_reactions:
            module_reactions[module_id] = [reaction_id]
        else:
            module_reactions[module_id].append(reaction_id)

    response_text = KEGG_BIOSERVICES.list('module')
    if response_text != '':
        csvreader = csv.reader(response_text.splitlines(), delimiter='\t')
    else:
        csvreader = []

    modules = {}
    for line in csvreader:
        module_id = line[0]
        module_data = line[1].split(', ')
        if len(module_data) > 2:
            module_name = module_data[0]
            module_formula = module_data[-1]
        else:
            module_name = module_data[0]
            module_formula = ''
        modules[module_id] = (module_name, module_formula)

    all_modules = set(list(modules.keys()) + list(module_reactions.keys()))
    with open(module_file, 'w') as output_file:
        csvwriter = csv.writer(output_file, delimiter='\t')
        csvwriter.writerow(['module_id', 'module_name', 'module_formula', 'module_reactions'])
        for module_id in all_modules:
            if module_id in module_reactions:
                module_reaction = ','.join(module_reactions[module_id])
            else:
                module_reaction = ''
            if module_id in modules:
                module_name = modules[module_id][0]
                module_formula = modules[module_id][1]
            else:
                module_name = ''
                module_formula = ''
            csvwriter.writerow([module_id, module_name, module_formula, module_reaction])


def create_sbml_model_from_kegg_file(reaction_folder, compound_file, output_sbml, output_tsv, pathways_tsv, remove_ubiquitous=True):
    """Using the reaction keg files (from retrieve_reactions), the compound file (from get_compound_names),
    create a SBML file (containing all reactions of KEGG) and a tsv file (used to map EC and/or KO to reaction ID)

    Args:
        reaction_folder (str): path to the folder containing keg reaction files
        compound_file (str): path to the tsv file containing compound ID and name
        output_sbml (str): path to the sbml output file
        output_tsv (str): path to an output tsv file mapping reaction ID, with KO and EC
        pathways_tsv (str): path to an output tsv showing the pathway, pathway ID and the associated reactions
    """
    compounds = {}
    with open(compound_file, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compounds[line[0]] = line[1]

    reaction_ecs = {}
    model = Model('KEGG')
    sbml_reactions = []
    reactions = {}
    pathways = {}
    for reaction_file in os.listdir(reaction_folder):
        reaction_file_path = os.path.join(reaction_folder, reaction_file)
        with open(reaction_file_path) as open_reaction_file_path:
            reaction_data = KEGG_BIOSERVICES.parse(open_reaction_file_path.read())

        reaction_id = reaction_data['ENTRY'].split(' ')[0]
        if 'NAME' in reaction_data:
            reaction_name = reaction_data['NAME'][0]
        else:
            reaction_name = None
        left_compounds, right_compounds = extract_reaction(reaction_id, reaction_data['EQUATION'])
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
        for stochiometry_metabolite in left_compounds:
            metabolite_id = stochiometry_metabolite[0]
            if remove_ubiquitous is True:
                if metabolite_id in UBIQUITOUS_METABOLITES:
                    continue
            metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=compounds[metabolite_id])
            reactions[reaction_id][metabolite_id_sbml] = - stochiometry_metabolite[1]
        for stochiometry_metabolite in right_compounds:
            metabolite_id = stochiometry_metabolite[0]
            if remove_ubiquitous is True:
                if metabolite_id in UBIQUITOUS_METABOLITES:
                    continue
            metabolite_id_sbml = Metabolite(metabolite_id, compartment='c', name=compounds[metabolite_id])
            reactions[reaction_id][metabolite_id_sbml] = stochiometry_metabolite[1]

        reaction = Reaction(reaction_id)
        if reaction_name:
            reaction.name = reaction_name
        if kegg_orthologs != []:
            reaction.gene_reaction_rule = '( ' + ' or '.join([ko for ko in kegg_orthologs]) + ' )'
        reaction.add_metabolites(reactions[reaction_id])

        sbml_reactions.append(reaction)

    model.add_reactions(sbml_reactions)
    logger.info('|kegg2bipartitegraph|reference| {0} reactions and {1} metabolites in reference model.'.format(len(model.reactions), len(model.metabolites)))
    # Create sbml file.
    write_sbml_model(model, output_sbml)

    with open(output_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['reaction_id', 'kegg_orthologs', 'kegg_ec'])
        for reaction_id in reaction_ecs:
            csvwriter.writerow([reaction_id, ','.join(reaction_ecs[reaction_id][0]), ','.join(reaction_ecs[reaction_id][1])])

    with open(pathways_tsv, 'w') as open_output_tsv:
        csvwriter = csv.writer(open_output_tsv, delimiter='\t')
        csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_reactions'])
        for pathway_id in pathways:
            pathway_name = pathways[pathway_id]['name']
            pathway_reactions = ','.join(set(pathways[pathway_id]['reaction']))
            csvwriter.writerow([pathway_id, pathway_name, pathway_reactions])


def get_kegg_database_version():
    """Retrieve the KEGG version release

    Returns:
        kegg_version (str): string containing KEGG release version
    """
    response_text = KEGG_BIOSERVICES.dbinfo()

    for line in response_text.splitlines():
        if 'Release ' in line:
            kegg_version = line.split('Release ')[1].strip()

    return kegg_version


def create_reference_base():
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|reference| Begin KEGG metabolism mapping.')

    # Create a json file containing metadata.
    options = {}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['bioservices'] = bioservices_version
    options['tool_dependencies']['python_package']['urllib'] = urllib.request.__version__
    options['tool_dependencies']['python_package']['cobra'] = cobra_version

    kegg2bipartitegraph_reference_metadata = {}
    kegg2bipartitegraph_reference_metadata['tool_options'] = options
    kegg2bipartitegraph_reference_metadata['kegg_release_number'] = get_kegg_database_version()

    output_folder = DATA_ROOT
    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(output_folder, 'kegg_model')
    is_valid_dir(kegg_model_path)

    kegg_reactions_folder_path = os.path.join(kegg_model_path, 'reaction_folder')
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')

    logger.info('|kegg2bipartitegraph|reference| Check missing files in {0}.'.format(DATA_ROOT))
    input_files = [kegg_sbml_model_path, kegg_rxn_mapping_path, kegg_pathways_path, kegg_modules_path]
    missing_files = []
    for input_file in input_files:
        if not os.path.exists(input_file):
            missing_files.append(input_file)
    
    if len(missing_files) > 0:
        logger.info('|kegg2bipartitegraph|reference| Missing: ' + ' '.join(missing_files))
        if not os.path.exists(kegg_reactions_folder_path):
            logger.info('|kegg2bipartitegraph|reference| Retrieve reactions from KEGG to create SMBL model.')
            get_reactions(kegg_reactions_folder_path)
        if not os.path.exists(compound_file_path):
            logger.info('|kegg2bipartitegraph|reference| Retrieve compound IDs and names from KEGG to create SMBL model.')
            get_compound_names(compound_file_path)
        logger.info('|kegg2bipartitegraph|reference| Create KEGG reference SBML and mapping tsv file.')
        create_sbml_model_from_kegg_file(kegg_reactions_folder_path, compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path, kegg_pathways_path)
        get_modules(kegg_modules_path)

        # Create compress archive.
        model_zipfile = zipfile.ZipFile(KEGG_ARCHIVE, mode="w")
        for filepath in input_files:
            file_name = os.path.basename(filepath)
            model_zipfile.write(filepath, file_name)

        model_zipfile.write(compound_file_path, 'kegg_compound_name.tsv')
        model_zipfile.close()
    else:
        logger.info('|kegg2bipartitegraph|reference| No missing files.')

    endtime = time.time()

    duration = endtime - starttime

    logger.info('|kegg2bipartitegraph|reference| Reference creation finished in {0}.'.format(duration))
