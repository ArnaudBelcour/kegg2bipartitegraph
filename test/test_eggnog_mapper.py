import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.eggnog import create_eggnog_network

def test_create_eggnog_network():
    create_eggnog_network('eggnog-mapper', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.tsv.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209', 'R03789', 'R01099', 'R01098']
    expected_metabolites = ['C00242', 'C01449', 'C01977', 'C01978', 'C16675', 'C20446', 'C00124', 'C00007', 'C03269', 'C00027', 'C00880']
    expected_modules_pathways_ids = ['map01100', 'map00052']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')

def test_create_eggnog_network_cli():
    subprocess.call(['k2bg', 'reconstruct_from_eggnog', '-i', 'eggnog-mapper', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.tsv.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209', 'R03789', 'R01099', 'R01098']
    expected_metabolites = ['C00242', 'C01449', 'C01977', 'C01978', 'C16675', 'C20446', 'C00124', 'C00007', 'C03269', 'C00027', 'C00880']
    expected_modules_pathways_ids = ['map01100', 'map00052']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')
