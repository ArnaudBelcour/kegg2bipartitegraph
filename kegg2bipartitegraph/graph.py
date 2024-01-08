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

import networkx as nx

import libsbml

def sbml_to_graphml(sbml_graph_file, graphml_file_path):
    """From a sbml file create a graphml file containing a bipartite graph.

    Args:
        sbml_graph_file (str): path to the sbml file
        graphml_file_path (str): path to the graphml
    """
    metabolic_graph = nx.DiGraph()

    reader = libsbml.SBMLReader()
    kegg_document = reader.readSBML(sbml_graph_file)
    kegg_model = kegg_document.getModel()

    for reaction in kegg_model.getListOfReactions():

        reactants = reaction.getListOfReactants()
        for reactant in reactants:
            metabolic_graph.add_edge(reactant.species, reaction.id)

        products = reaction.getListOfProducts()
        for product in products:
            metabolic_graph.add_edge(reaction.id, product.species)

    nx.write_graphml(metabolic_graph, graphml_file_path)