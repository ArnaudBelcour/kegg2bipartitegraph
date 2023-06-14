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

from kegg2bipartitegraph.utils import is_valid_dir
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.graph import sbml_to_graphml

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)

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
    if '<=>' not in equation_text:
        logger.critical('|kegg2bipartitegraph|reference| No <=> symbol in equation {0} of {1}.'.format(equation_text, reaction_id))
    if len(equation_text.split('<=>')) > 2:
        logger.critical('|kegg2bipartitegraph|reference| More than one <=> symbol in equation {0} of {1}.'.format(equation_text, reaction_id))

    # Find group pattern
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
                logger.critical('|kegg2bipartitegraph|reference| Stochiometry is not an int ({0}) for {1}, replace by 1 as it is not relevant for topological analysis.'.format(stoechiometry, reaction_id))
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


def libsbml_check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value == None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + libsbml.OperationReturnValue_toString(value).strip() + '"'
            raise TypeError(err_msg)
    else:
        return


def create_sbml_model_from_kegg_file_libsbml(reaction_folder, compound_file, output_sbml, output_tsv, pathways_tsv,
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
        pathways_tsv (str): path to an output tsv showing the pathway, pathway ID and the associated reactions
        remove_ubiquitous (bool): remove ubiquitous metabolites
        remove_glycan_reactions (bool): remove glycan associated reactions
    """
    # Set libsbml namespace.
    sbml_ns = libsbml.SBMLNamespaces(3, 1)  # SBML L3V1
    sbml_ns.addPackageNamespace("fbc", 2)  # fbc-v2

    # Create document
    document = libsbml.SBMLDocument(sbml_ns)
    model = document.createModel('KEGG')

    document.enablePackage(libsbml.FbcExtension.getXmlnsL3V1V2(), 'fbc', True)
    document.setPackageRequired("fbc", False)
    model_fbc = model.getPlugin('fbc')
    model_fbc.setStrict(True)

    # Set units.
    libsbml_check(model,                              'create model')
    libsbml_check(model.setTimeUnits("second"),       'set model-wide time units')
    libsbml_check(model.setExtentUnits("mole"),       'set model units of extent')
    libsbml_check(model.setSubstanceUnits('mole'),    'set model substance units')

    math_ast = libsbml.parseL3Formula('FLUX_VALUE')
    libsbml_check(math_ast, 'create AST for rate expression')

    # Set compartments.
    compart = model.createCompartment()
    libsbml_check(compart,'create compartment')
    libsbml_check(compart.setId('c'),'set compartment id c')
    libsbml_check(compart.setSize(1),'set size for compartment id c')
    libsbml_check(compart.setConstant(True),'set constant for compartment id c')
    libsbml_check(compart.setName("cytosol"),'set compartment name cytosol')

    # Set default bound values.
    default_lb = model.createParameter()
    default_lb.setId('default_lower_bound')
    default_lb.setValue(-1000)
    default_lb.setConstant(True)

    default_ub = model.createParameter()
    default_ub.setId('default_upper_bound')
    default_ub.setValue(1000)
    default_ub.setConstant(True)

    zero_bound = model.createParameter()
    zero_bound.setId('default_zero_bound')
    zero_bound.setValue(0)
    zero_bound.setConstant(True)

    # Read compounds file.
    compounds = {}
    with open(compound_file, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compounds[line[0]] = line[1]

    reaction_ecs = {}
    sbml_reactions = []
    reactions = {}
    genes = []
    pathways = {}
    metabolites = []

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
        for stochiometry_metabolite in left_compounds:
            metabolite_id = stochiometry_metabolite[0]
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
            libsbml_check(species_ref.setStoichiometry(stochiometry_metabolite[1]), 'set stoichiometry {0}'.format(stochiometry_metabolite[1]))
            libsbml_check(species_ref.setConstant(False), 'set constant %s' %False)
            reactants.append(metabolite_id)

        products = []
        # Create metabolites from right compounds, remove ubiquitous and mark glycan reactions.
        for stochiometry_metabolite in right_compounds:
            metabolite_id = stochiometry_metabolite[0]
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
            libsbml_check(species_ref.setStoichiometry(stochiometry_metabolite[1]), 'set stoichiometry {0}'.format(stochiometry_metabolite[1]))
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

        remove_reaction = False
        reactions_metabolites = reactants + products
        if reactions_metabolites == []:
            logger.critical('|kegg2bipartitegraph|reference| No reactants and products for {0}, will be removed from model.'.format(reaction_id))
            remove_reaction = True
        if remove_glycan_reactions is True and glycan_reaction is True:
            logger.critical('|kegg2bipartitegraph|reference| Do not add glycan reaction {0}.'.format(reaction_id))
            remove_reaction = True

        # Remove reaction if it contains glycan metabolite or it does not have reactants and products.
        if remove_reaction is True:
            model.removeReaction(reaction_id)

    libsbml.writeSBMLToFile(document, output_sbml)
    logger.info('|kegg2bipartitegraph|reference| {0} reactions and {1} metabolites in reference model.'.format(len(model.getListOfReactions()), len(model.getListOfSpecies())))

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
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__

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
    kegg_graphml_model_path = os.path.join(kegg_model_path, 'kegg_model.graphml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_metadata_path = os.path.join(kegg_model_path, 'kegg_metadata.json')

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
        create_sbml_model_from_kegg_file_libsbml(kegg_reactions_folder_path, compound_file_path, kegg_sbml_model_path, kegg_rxn_mapping_path, kegg_pathways_path)
        sbml_to_graphml(kegg_sbml_model_path, kegg_graphml_model_path)
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

    kegg2bipartitegraph_reference_metadata['kegg2bipartitegraph_reference_duration'] = duration
    with open(kegg_metadata_path, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_reference_metadata, ouput_file, indent=4)

    logger.info('|kegg2bipartitegraph|reference| Reference creation finished in {0}.'.format(duration))
