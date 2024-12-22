import shutil
import subprocess
import libsbml
import os
import json

from kegg2bipartitegraph.scope import compute_scope

def test_compute_scope():
    graphml_path = os.path.join('scope_data', 'graphml')
    seed_path = os.path.join('scope_data', 'scope_seed.txt')
    output_folder = 'test_out'
    compute_scope(graphml_path, output_folder, seed_path)

    output_json_file = os.path.join(output_folder, 'accessibility.json')
    with open(output_json_file, 'r') as open_output_json_file:
        json_data = json.load(open_output_json_file)

    expected_scope = ['C_A', 'C_B', 'C_C', 'C_D', 'C_E', 'C_F', 'C_G', 'C_H']
    expected_activated_reactions = ['R1', 'R2', 'R3']

    found_scope = json_data['producible_metabolites']['orgA']
    found_activated_reactions = json_data['activated_reactions']['orgA']

    assert sorted(expected_scope) == sorted(found_scope)
    assert sorted(expected_activated_reactions) == sorted(found_activated_reactions)

    shutil.rmtree('test_out')

def test_compute_scope_cli():
    graphml_path = os.path.join('scope_data', 'graphml')
    seed_path = os.path.join('scope_data', 'scope_seed.txt')
    output_folder = 'test_out'
    subprocess.call(['k2bg', 'scope', '-i', graphml_path, '-o', output_folder, '-s', seed_path])

    output_json_file = os.path.join(output_folder, 'scope_seed', 'accessibility.json')
    with open(output_json_file, 'r') as open_output_json_file:
        json_data = json.load(open_output_json_file)

    expected_scope = ['C_A', 'C_B', 'C_C', 'C_D', 'C_E', 'C_F', 'C_G', 'C_H']
    expected_activated_reactions = ['R1', 'R2', 'R3']

    found_scope = json_data['producible_metabolites']['orgA']
    found_activated_reactions = json_data['activated_reactions']['orgA']

    assert sorted(expected_scope) == sorted(found_scope)
    assert sorted(expected_activated_reactions) == sorted(found_activated_reactions)

    shutil.rmtree('test_out')
