import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.kofamkoala import create_kofamkoala_network

def test_create_kofamkoala_network():
    create_kofamkoala_network('kofam_koala', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209', 'R01098']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242', 'C00124', 'C00007', 'C00027', 'C00880']
    expected_modules_pathways_ids = ['map01100', 'map00052']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')

def test_create_kofamkoala_network_cli():
    subprocess.call(['k2bg', 'reconstruct_from_kofamkoala', '-i', 'kofam_koala', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209', 'R01098']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242', 'C00124', 'C00007', 'C00027', 'C00880']
    expected_modules_pathways_ids = ['map01100', 'map00052']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')