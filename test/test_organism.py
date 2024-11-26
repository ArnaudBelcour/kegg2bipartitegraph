import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.organism import create_organism_network

def test_create_organism_network():
    create_organism_network('eco', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    expected_len_reactions = 1771
    expected_len_metabolites = 1760
    expected_len_pathways = 100
    expected_len_modules = 181

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    pathways_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('map')]
    modules_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('M')]

    assert len(found_reactions) == expected_len_reactions
    assert len(found_metabolites) == expected_len_metabolites
    assert len(pathways_ids) == expected_len_pathways
    assert len(modules_ids) == expected_len_modules

    shutil.rmtree('test_out')

def test_create_organism_network_cli():
    subprocess.call(['k2bg', 'reconstruct_from_organism', '-i', 'eco', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    expected_len_reactions = 1771
    expected_len_metabolites = 1760
    expected_len_pathways = 100
    expected_len_modules = 181

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    pathways_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('map')]
    modules_ids = [group.id for group in model_groups.getListOfGroups() if group.id.startswith('M')]

    assert len(found_reactions) == expected_len_reactions
    assert len(found_metabolites) == expected_len_metabolites
    assert len(pathways_ids) == expected_len_pathways
    assert len(modules_ids) == expected_len_modules

    shutil.rmtree('test_out')