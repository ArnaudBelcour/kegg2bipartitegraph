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

import csv
import logging
import os
import urllib.parse
import urllib.request
import libsbml

from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version

URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)


def get_go_to_ec(ec_to_gos_file):
    """From the EC to GO tsv mapping file (creating at go_to_ec in reference file)
    retrieve a mapping dictionary between EC and GO

    Args:
        ec_to_gos_file (str): path to the EC to GO tsv file

    Returns:
        go_to_ecs (dict): mapping between GO term and list of EC numbers as value
    """
    go_to_ecs = {}
    with open(ec_to_gos_file, 'r') as open_ec_to_gos_file:
        csvreader = csv.reader(open_ec_to_gos_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            ec_id = line[0]
            gos = line[1].split(',')
            for go in gos:
                if go not in go_to_ecs:
                    go_to_ecs[go] = [ec_id]
                else:
                    go_to_ecs[go].append(ec_id)
    return go_to_ecs


def retrieve_mapping_dictonaries(kegg_rxn_mapping_path):
    """From the KEGG tsv mapping file (creating at create_sbml_model_from_kegg_file)
    retrieve 2 mapping dictionaries, one KO ID to kegg_reaction and another
    EC_number to kegg_reaction

    Args:
        kegg_rxn_mapping_path (str): path to the KEGG tsv file for mapping reaction and EC/KO

    Returns:
        ko_to_reactions (dict): mapping between KO ID as key and kegg reaction as value
        ec_to_reactions (dict): mapping between EC number as key and kegg reaction as value
    """
    ko_to_reactions = {}
    ec_to_reactions = {}
    with open(kegg_rxn_mapping_path, 'r') as input_file:
        csvreader = csv.reader(input_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            reaction_id = line[0]
            if line[1] != '':
                ko_ids = line[1].split(',')
            else:
                ko_ids = []
            if line[2] != '':
                ec_ids = line[2].split(',')
            else:
                ec_ids = []
            for ko_id in ko_ids:
                if ko_id not in ko_to_reactions:
                    ko_to_reactions[ko_id] = [reaction_id]
                else:
                    ko_to_reactions[ko_id].append(reaction_id)

            for ec_id in ec_ids:
                if ec_id not in ec_to_reactions:
                    ec_to_reactions[ec_id] = [reaction_id]
                else:
                    ec_to_reactions[ec_id].append(reaction_id)

    return ko_to_reactions, ec_to_reactions


def compute_stat_kegg(sbml_folder, stat_file=None):
    """Compute stat associated with the number of reactions and metabolites for each taxonomic affiliations.

    Args:
        sbml_folder (str): pathname to the sbml folder containing draft metabolic networks for each cluster
        stat_file (str): pathname to the tsv stat file

    Returns:
        kegg_numbers (dict): dict containing observation names (as key) associated with reaction and metabolites number (as value)
    """
    kegg_numbers = {}
    for infile in os.listdir(sbml_folder):
        if '.sbml' in infile:
            sbml_input_file_path = os.path.join(sbml_folder, infile)
            reader = libsbml.SBMLReader()
            kegg_document = reader.readSBML(sbml_input_file_path)
            kegg_model = kegg_document.getModel()
            infile_reactions = kegg_model.getListOfReactions()
            infile_metabolites = kegg_model.getListOfSpecies()
            kegg_numbers[infile.replace('.sbml','')] = (len(infile_reactions), len(infile_metabolites))

    if stat_file:
        with open(stat_file, 'w') as stat_file_open:
            csvwriter = csv.writer(stat_file_open, delimiter='\t')
            csvwriter.writerow(['observation_name', 'Number_reactions', 'Number_metabolites'])
            for observation_name in kegg_numbers:
                csvwriter.writerow([observation_name, kegg_numbers[observation_name][0], kegg_numbers[observation_name][1]])

    return kegg_numbers