from kegg2bipartitegraph.reference import extract_reaction

def test_extract_reaction_simple():
    expected_left_compounds = [('C00479', 2)]
    expected_right_compounds = [('C02948', 1)]
    reaction_id = 'R00038'
    equation_text = '2 C00479 <=> C02948'
    left_compounds, right_compounds, modify_stoichiometry= extract_reaction(reaction_id, equation_text)

    assert expected_left_compounds == left_compounds
    assert expected_right_compounds == right_compounds


def test_extract_reaction_n():
    expected_left_compounds = [('C02321', 1), ('C00007', 1)]
    expected_right_compounds = [('C21834', 1)]
    reaction_id = 'R11999'
    equation_text = 'C02321 + n C00007 <=> n C21834'

    left_compounds, right_compounds, modify_stoichiometry = extract_reaction(reaction_id, equation_text)

    assert expected_left_compounds == left_compounds
    assert expected_right_compounds == right_compounds