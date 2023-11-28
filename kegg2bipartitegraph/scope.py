# Copyright (C) 2021-2023 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


def test_seeds_graph(seeds, model_graph) :
    """ Tests if seeds exist in graph, else only selects compounds in graph,
    otherwise throws an error """
    if not np.all(np.in1d(list(seeds), list(model_graph.nodes)) == True):
        length = len(seeds)
        logger.warning("At least one wrong compound name for inputs, will be removed")
        seeds = list(np.array(seeds)[np.in1d(list(seeds), list(model_graph.nodes))])
        logger.warning("%d/%d kept from added inputs" %(len(seeds), length))
        if len(seeds) < 1 :
            logger.error("Not enough inputs")
            raise SystemExit()
    return seeds


def compute_scope_from_graphml(graphml_file, seeds):
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
