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

import argparse
import csv
import datetime
import json
import logging
import os
import urllib.request
import sys
import time

from socket import timeout

from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version

MIN_VAL = 0
MAX_VAL = 1
URLLIB_HEADERS = {'User-Agent': 'kegg2bipartitegraph annotation v' + kegg2bipartitegraph_version + ', request by urllib package v' + urllib.request.__version__}

logger = logging.getLogger(__name__)


def range_limited_float_type(arg):
    """Type function for argparse - a float within some predefined bounds

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")
    if f < MIN_VAL or f > MAX_VAL:
        raise argparse.ArgumentTypeError("Argument must be < " + str(MAX_VAL) + " and > " + str(MIN_VAL))
    return f


def limited_integer_type(arg):
    """Type function for argparse - an integer

    Args:
        arg: argparse argument

    Returns:
        arg: argparse argument
    """
    try:
        f = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be an integer number")
    return f


def is_valid_path(filepath):
    """Return True if filepath is valid
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it)
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True


def urllib_query(request, nb_retry=5):
    """Use urllib to query UniProt.

    Args:
        request (urllib.request.Request): request object
        nb_retry (int): number of retry to perform (default = 5)

    Returns:
        response (http.client.HTTPResponse): response returns by urllib.request
    """
    passed = False

    if nb_retry == 0:
        sys.exit('5 retry attempts have been performed but were not successful, so esmecata has been stopped. You may try to relaunch it using the same command esmecata should resume.')
    try:
        response = urllib.request.urlopen(request)
        passed = True
    except timeout:
        logger.critical('Timeout occurs for query, try to relaunch query.')
        time.sleep(10)
        urllib_query(request, nb_retry-1)

    if passed is True:
        return response


def get_rest_uniprot_release(options):
    """Get the release version and date of Uniprot and Trembl and also the date of the query.

    Args:
        options (dict): dictionary with metadata on the esmecata run

    Returns:
        dict: metadata of Uniprot release
    """
    uniprot_releases = {}

    # Get Uniprot release version
    uniprot_urllib_request = urllib.request.Request('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt',
                                                    headers=URLLIB_HEADERS)
    uniprot_response = urllib_query(uniprot_urllib_request)
    uniprot_lines = uniprot_response.readlines()
    uniprot_release_number = uniprot_lines[0].decode('utf-8').split(' ')[3].replace('\n','')
    swissprot_release_number = uniprot_lines[1].decode('utf-8').split(' ')[2].replace('\n','')
    swissprot_release_date = uniprot_lines[1].decode('utf-8').split(' ')[4].replace('\n','')
    trembl_release_number = uniprot_lines[2].decode('utf-8').split(' ')[2].replace('\n','')
    trembl_release_date = uniprot_lines[2].decode('utf-8').split(' ')[4].replace('\n','')

    date = datetime.datetime.now().strftime('%d-%B-%Y %H:%M:%S')

    uniprot_releases['esmecata_query_system'] = 'REST queries on Uniprot'
    uniprot_releases['uniprot_release'] = uniprot_release_number
    uniprot_releases['access_time'] = date
    uniprot_releases['swissprot_release_number'] = swissprot_release_number
    uniprot_releases['swissprot_release_date'] = swissprot_release_date
    uniprot_releases['trembl_release_number'] = trembl_release_number
    uniprot_releases['trembl_release_date'] = trembl_release_date
    uniprot_releases['tool_options'] = options

    return uniprot_releases


def get_internal_reference_model_path():
    """ Retrieve path to reference model.

    Returns:
        kegg_model_path (str): Path to reference model
    """
    reference_root = os.path.dirname(__file__)
    kegg_model_path = os.path.join(reference_root, 'data', 'kegg_model')

    return kegg_model_path


def get_current_database_version(database_path):
    """ Retrieve database version, using metadata file.

    Returns:
        database_version (str): KEGG database version of reference model
    """
    kegg_json_model_path = os.path.join(database_path, 'kegg_metadata.json')

    with open(kegg_json_model_path, 'r') as input_metadata_json:
        json_data = json.load(input_metadata_json)
    database_version = json_data['kegg_release_number']

    return database_version


def write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions):
    """ From KEGG pathways and reactions found in an organism, write a pathway file.

    Args:
        kegg_pathways (dict): dictionary with pathway ID as key and a tuple with pathway name and pathway reactions as value
        pathways_output_file_path (str): path to pathway output file
        total_added_reactions (list): list of reactions found in organism
    """
    organism_pathways = []

    with open(pathways_output_file_path, 'w') as open_pathways_output_file_path:
        csvwriter = csv.writer(open_pathways_output_file_path, delimiter='\t')
        csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_completion_ratio', 'pathway_reaction_in_taxon', 'pathway_reaction'])
        for pathway in kegg_pathways:
            pathway_reactions = kegg_pathways[pathway][1]
            pathway_reaction_in_taxon = set(pathway_reactions).intersection(set(total_added_reactions))
            if len(pathway_reaction_in_taxon) > 0:
                pathway_name = kegg_pathways[pathway][0]
                pathway_completion_ratio = len(pathway_reaction_in_taxon) / len(pathway_reactions)
                if pathway_completion_ratio > 0:
                    organism_pathways.append(pathway)
                    csvwriter.writerow([pathway, pathway_name, pathway_completion_ratio, ','.join(pathway_reaction_in_taxon), ','.join(pathway_reactions)])

    return organism_pathways


def write_module_file(kegg_modules, modules_output_file_path, total_added_reactions, min_nb_rxn=0, min_ratio_rxn=0):
    """ From KEGG modules and reactions found in an organism, write a pathway file.

    Args:
        kegg_modules (dict): dictionary with module ID as key and a tuple with module name and module reactions as value
        modules_output_file_path (str): path to module output file
        total_added_reactions (list): list of reactions found in organism
        min_nb_rxn (int): minimal number of reactions to keep a pathway
        min_ratio_rxn (float): minimal ratio of reactions to keep a pathway
    """
    organism_modules = []

    with open(modules_output_file_path, 'w') as open_modules_output_file_path:
        csvwriter = csv.writer(open_modules_output_file_path, delimiter='\t')
        csvwriter.writerow(['module_id', 'module_name', 'module_completion_ratio', 'module_reaction_in_taxon', 'module_reaction'])
        for module in kegg_modules:
            module_reactions = kegg_modules[module][1]
            module_reaction_in_taxon = set(module_reactions).intersection(set(total_added_reactions))
            if len(module_reaction_in_taxon) > min_nb_rxn:
                module_name = kegg_modules[module][0]
                module_completion_ratio = len(module_reaction_in_taxon) / len(module_reactions)
                if module_completion_ratio > min_ratio_rxn:
                    organism_modules.append(module)
                    csvwriter.writerow([module, module_name, module_completion_ratio, ','.join(module_reaction_in_taxon), ','.join(module_reactions)])

    return organism_modules