# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

# This code is an adaptation of the scope developed by Adèle Webber Zendrera in:
# https://doi.org/10.1038/s41598-021-91486-8
# It is available here:
# https://github.com/AWebZen/FunctionalPrediction5000species/blob/2ff7c4fd4092a8b565da9d3fa2f4b557d9954bf0/utils_objet.py#L99-L170
# As the purpose of this package is to generalise this approach to be able to reconstruct metabolic networks, I have also added the scope method.
# WARNING: I have modified the method so it does not work with matrix and graphml associated with the source article.

import logging
import time
import os
import csv
import json
import sys
import networkx as nx
import numpy as np

from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.utils import is_valid_dir

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')


def test_seeds_graph(seeds, model_graph) :
    """Function from the work of Adèle Webber Zendrera et al. (2021),
    Tests if seeds exist in graph, else only selects compounds in graph,
    otherwise throws an error.

    Args:
        seeds (list): list of seeds (KEGG IDs)
        model_graph (networkx object): networkx object containing metabolic bipartite graph

    Returns:
        seeds (list): list of seeds in the graph
    """
    if not np.all(np.in1d(list(seeds), list(model_graph.nodes)) == True):
        length = len(seeds)
        seeds = list(np.array(seeds)[np.in1d(list(seeds), list(model_graph.nodes))])

        logger.warning("|kegg2bipartitegraph|scope| At least one metabolite absent from graph (or wrong compound name), will be removed: {0}/{1} kept from added inputs.".format(len(seeds), length))
        if len(seeds) < 1 :
            logger.error("|kegg2bipartitegraph|scope| Not enough inputs")
            raise SystemExit()
    return seeds


def compute_scope_from_graphml(graphml_file, seeds):
    """Function modified from the work of Adèle Webber Zendrera et al. (2021),
    compute the scope using the list of seeds and metaoblic bipartite graph.

    Args:
        graphml_file (str): path to the input graphml file
        seeds (list): list of seeds (KEGG IDs)

    Returns:
        accessibility (dict): dictionary showing the accessible and not accessible reactions and metabolites
    """
    # Read the bipartite graph in graphml format. 
    model_graph = nx.read_graphml(graphml_file)
    # Remove isolated nodes
    model_graph.remove_nodes_from(list(nx.isolates(model_graph)))

    # Source compounds (seeds) from which accessibility will be measured
    if seeds:
        seeds = test_seeds_graph(seeds, model_graph)
    else:
        # Inputs deduced from graph ("complete" input compounds)
        in_array = np.array(list(model_graph.in_degree)) #reaction graph's in-degrees
        input_cpd = in_array[np.where(in_array[:,1] == '0')[0],0] #out nodes of in/out graph (sources)
        if added_inp :
            added_inp = test_seeds_graph(added_inp)
        seeds = set(list(input_cpd) + added_inp) #+ ["C00028", "C00073"]

    accessibility = dict.fromkeys(model_graph.nodes, "Non accessible")

    #Initialisation: all seeds are accessible
    for seed in seeds:
        accessibility[seed] = "Accessible"
        if seed.startswith('R'):
            # If the seed node is a reaction node, make accessible all of its sucessords (= products).
            successors = list(model_graph.successors(seed))
            for successor in successors:
                accessibility[successor] = "Accessible"

    #Accessible nodes from each seed.
    for seed in seeds:
        for _, traversed_node in nx.bfs_edges(model_graph, seed):
            # Select only reaction node.
            if not traversed_node.startswith('R'):
                continue
            # Check if all predecessors of reaction node (= reactants) are accessible.
            predecessors = list(model_graph.predecessors(traversed_node))
            acessible_predecessors = np.array([accessibility[predecessor] for predecessor in predecessors])
            if np.all(acessible_predecessors == "Accessible"):
                # If yes, all sucessors of reaction node (= products) are accessible.
                accessibility[traversed_node] = "Accessible"
                successors = list(model_graph.successors(traversed_node))
                for successor in successors :
                    accessibility[successor] = "Accessible"

    return accessibility


def compute_scope(input_graphml_folder, output_folder, seed_file, reference_folder=False):
    """Using graphml files created by k2bg compute scope using seed file.

    Args:
        input_graphml_folder (str): path to the input folder containing graphml files
        output_folder (str): path to the output folder
        seed_file (str): path to txt file containing seeds
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|scope| Begin scope computation for grpahml filesi n {0}.'.format(input_graphml_folder))

    if reference_folder is not False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|eggnog| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|eggnog| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Download Uniprot metadata and create a json file containing them.
    options = {'input_graphml_folder': input_graphml_folder, 'output_folder': output_folder, 'seed_file': seed_file}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['networkx'] = nx.__version__

    kegg2bipartitegraph_scope_metadata = {}
    kegg2bipartitegraph_scope_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    # Get compound names.
    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    compound_names = {}
    with open(compound_file_path, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compound_names[line[0]] = line[1]

    seed_metabolites = []
    with open(seed_file, 'r') as open_seed_file:
        for line in open_seed_file.readlines():
            seed_metabolites.append(line.strip('\n'))

    producible_metabolites = {}
    activated_reactions = {}
    for graphml_file in os.listdir(input_graphml_folder):
        organism_name = graphml_file.replace('.graphml', '')
        logger.info('|kegg2bipartitegraph|scope| -- Check producibility of {0}'.format(organism_name))
        graphml_path = os.path.join(input_graphml_folder, graphml_file)
        accessibility = compute_scope_from_graphml(graphml_path, seed_metabolites)
        producible_compounds = [compound_names[compound] for compound in accessibility if accessibility[compound] == 'Accessible' and compound.startswith('C') and compound not in seed_metabolites]
        reachable_reactions = [reaction for reaction in accessibility if accessibility[reaction] == 'Accessible' and reaction.startswith('R')]
        producible_metabolites[organism_name] = producible_compounds
        activated_reactions[organism_name] = reachable_reactions
        logger.info('|kegg2bipartitegraph|scope| {0} producible metabolites.'.format(len(producible_compounds)))

    result_json = {}
    result_json['producible_metabolites'] = producible_metabolites
    result_json['activated_reactions'] = activated_reactions

    output_json_file = os.path.join(output_folder, 'accessibility.json')
    with open(output_json_file, 'w') as open_output_json_file:
        json.dump(result_json, open_output_json_file, indent=4)

    endtime = time.time()

    duration = endtime - starttime
    kegg2bipartitegraph_scope_metadata['kegg2bipartitegraph_esmecata_duration'] = duration
    kegg2bipartitgraph_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_esmecata_kegg.json')
    with open(kegg2bipartitgraph_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_scope_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|scope| Computation of scope analysis of metabolic networks done.')
