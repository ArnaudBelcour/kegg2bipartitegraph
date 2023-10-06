import shutil
import subprocess

from cobra.io import read_sbml_model

from kegg2bipartitegraph.esmecata import create_esmecata_network

def test_draft_reconstruct():
    create_esmecata_network('esmecata_annotation_folder', 'test_out')

    sbml_model = read_sbml_model('test_out/sbml/taxon_test.sbml')

    expected_reactions = ['R10209']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242']

    found_reactions = [reaction.id for reaction in sbml_model.reactions]
    found_metabolites = [reaction.id for reaction in sbml_model.metabolites]
    assert found_reactions == expected_reactions
    assert sorted(found_metabolites) == sorted(expected_metabolites)

    shutil.rmtree('test_out')

def test_draft_reconstruct_cli():
    subprocess.call(['kegg2bipartitegraph', 'reconstruct_from_esmecata', '-i', 'esmecata_annotation_folder', '-o', 'test_out'])
    sbml_model = read_sbml_model('test_out/sbml/taxon_test.sbml')

    expected_reactions = ['R10209']
    expected_metabolites = ['C01977', 'C16675', 'C20446', 'C00242']

    found_reactions = [reaction.id for reaction in sbml_model.reactions]
    found_metabolites = [reaction.id for reaction in sbml_model.metabolites]
    assert found_reactions == expected_reactions
    assert sorted(found_metabolites) == sorted(expected_metabolites)

    shutil.rmtree('test_out')