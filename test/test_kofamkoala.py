import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.kofamkoala import create_kofamkoala_network

def test_draft_reconstruct():
    create_kofamkoala_network('kofam_koala', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]
    assert found_reactions == expected_reactions
    assert sorted(found_metabolites) == sorted(expected_metabolites)

    shutil.rmtree('test_out')

def test_draft_reconstruct_cli():
    subprocess.call(['kegg2bipartitegraph', 'reconstruct_from_kofamkoala', '-i', 'kofam_koala', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/result.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R10209']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]
    assert found_reactions == expected_reactions
    assert sorted(found_metabolites) == sorted(expected_metabolites)

    shutil.rmtree('test_out')