from kegg2bipartitegraph.reference import extract_reaction, extract_parent_child_nested_dict

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


def test_extract_parent_child_nested_dict():
    input_dicitonary = {'A': {'A1': {'Aa1': ['element_a1', 'element_a2']}, 'A2': {'Aa2': ['element_aa1', 'element_a2']}},
                        'B': {'B2': {'Bb2': ['element_b2']}},
                        'C': {'C1': ['element_c1']}}

    parent_child_dict = {}
    extract_parent_child_nested_dict(input_dicitonary, parent_child_dict)
    expected_parent_child_dict = {'element_a1': ['Aa1'], 'element_a2': ['Aa1', 'Aa2'], 'element_aa1': ['Aa2'], 'element_b2': ['Bb2'], 'element_c1': ['C1'],
                                  'Aa1': ['A1'], 'Aa2': ['A2'], 'Bb2': ['B2'],
                                  'A1': ['A'], 'A2': ['A'], 'B2': ['B'], 'C1': ['C']}

    for child_element in parent_child_dict:
        assert sorted(parent_child_dict[child_element]) == sorted(expected_parent_child_dict[child_element])