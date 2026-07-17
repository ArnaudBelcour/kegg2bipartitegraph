import shutil
import subprocess
import libsbml
import csv
import os

import kegg2bipartitegraph
from kegg2bipartitegraph.organism import create_organism_network

kegg2bipartitegraph_path = kegg2bipartitegraph.__path__[0]
DATA_ROOT = os.path.join(kegg2bipartitegraph_path, 'data')
KEGG_MODEL_PATH = os.path.join(DATA_ROOT, 'kegg_model')
KEGG_REMOVED_CHANGED_REACTION_PATH = os.path.join(KEGG_MODEL_PATH, 'kegg_removed_changed_reaction.tsv')

EXPECTED_LEN_REACTIONS = 1825
EXPECTED_LEN_METABOLITES = 1765
EXPECTED_LEN_PATHWAYS = 105
EXPECTED_LEN_MODULES = 108

# Remove reactions associated with glycan.
reaction_to_removes = []
with open(KEGG_REMOVED_CHANGED_REACTION_PATH, 'r') as open_KEGG_REMOVED_CHANGED_REACTION_PATH:
    csvreader = csv.DictReader(open_KEGG_REMOVED_CHANGED_REACTION_PATH, delimiter='\t')
    for line in csvreader:
        reaction_to_removes.append(line['reaction_id'])

def test_create_organism_network():
    create_organism_network('eco', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    pathways_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('map')]
    modules_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('M')]

    # Check that no reactions to remove are present in the model.
    assert len(set(reaction_to_removes).intersection(set(found_reactions))) == 0
    assert len(found_reactions) == EXPECTED_LEN_REACTIONS
    assert len(found_metabolites) == EXPECTED_LEN_METABOLITES
    assert len(pathways_ids) == EXPECTED_LEN_PATHWAYS
    assert len(modules_ids) == EXPECTED_LEN_MODULES

    shutil.rmtree('test_out')

def test_create_organism_network_cli():
    subprocess.call(['k2bg', 'reconstruct_from_organism', '-i', 'eco', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    pathways_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('map')]
    modules_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('M')]

    # Check that no reactions to remove are present in the model.
    assert len(set(reaction_to_removes).intersection(set(found_reactions))) == 0
    assert len(found_reactions) == EXPECTED_LEN_REACTIONS
    assert len(found_metabolites) == EXPECTED_LEN_METABOLITES
    assert len(pathways_ids) == EXPECTED_LEN_PATHWAYS
    assert len(modules_ids) == EXPECTED_LEN_MODULES

    shutil.rmtree('test_out')