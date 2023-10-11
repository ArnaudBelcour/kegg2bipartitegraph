import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.organism import create_organism_network

def test_draft_reconstruct():
    create_organism_network('eco', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    expected_len_reactions = 1766
    expected_len_metabolites = 1756

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    assert len(found_reactions) == expected_len_reactions
    assert len(found_metabolites) == expected_len_metabolites

    shutil.rmtree('test_out')

def test_draft_reconstruct_cli():
    subprocess.call(['kegg2bipartitegraph', 'reconstruct_from_organism', '-i', 'eco', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/eco.sbml')
    sbml_model = sbml_document.getModel()

    expected_len_reactions = 1766
    expected_len_metabolites = 1756

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    assert len(found_reactions) == expected_len_reactions
    assert len(found_metabolites) == expected_len_metabolites

    shutil.rmtree('test_out')