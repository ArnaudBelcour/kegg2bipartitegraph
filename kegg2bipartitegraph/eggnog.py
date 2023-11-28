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

import csv
import logging
import os
import json
import time
import numpy as np
import itertools
import pandas as pd
import sys
import libsbml

from kegg2bipartitegraph.utils import is_valid_dir, write_pathway_file, write_module_file
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def read_annotation(eggnog_outfile:str, remove_duplicates=True):
    """Read an eggnog-mapper annotation file and retrieve EC numbers and GO terms by genes.

    Args:
        eggnog_outfile (str): path to eggnog-mapper annotation file
        remove_duplicates (str): if there are multiple time the same gene ID, take the one with the highest score
    Returns:
        dict: dict of genes and their annotations as {gene1:{EC:'..,..', GOs:'..,..,'}}
    """
    # Look at the twentieth first rows to find the header.
    with open(eggnog_outfile, 'r') as f:
        twentieth_first_rows = list(itertools.islice(f, 20))
        first_row_after_header = min([index for index, str_row in enumerate(twentieth_first_rows) if not str_row.startswith('#')])
        header_row = first_row_after_header - 1
        headers_row = twentieth_first_rows[header_row].lstrip("#").strip().split('\t')

    # Fix issue when header is incomplete (eggnog before version 2.0).
    if len(headers_row) == 17:
        headers_row.extend(['tax_scope', 'eggNOG_OGs', 'bestOG', 'COG_functional_category', 'eggNOG_free_text'])

    to_extract_annotations = ['GOs','EC', 'Preferred_name']
    if 'PFAMs' in headers_row:
        to_extract_annotations.append('PFAMs')
    if 'BiGG_Reaction' in headers_row:
        to_extract_annotations.append('BiGG_Reaction')
    if 'KEGG_Reaction' in headers_row:
        to_extract_annotations.append('KEGG_Reaction')
    if 'KEGG_ko' in headers_row:
        to_extract_annotations.append('KEGG_ko')
    if 'CAZy' in headers_row:
        to_extract_annotations.append('CAZy')

    # Use chunk when reading eggnog file to cope with big file.
    chunksize = 10 ** 6
    for annotation_data in pd.read_csv(eggnog_outfile, sep='\t', comment='#', header=None, dtype = str, chunksize = chunksize):
        annotation_data.replace(np.nan, '', inplace=True)
        # Assign the headers
        annotation_data.columns = headers_row

        if remove_duplicates is True:
            # If there are multiple times the same gene ID, select the one with the best score.
            annotation_data = annotation_data.sort_values('score', ascending=False).drop_duplicates('query').sort_index()
        if 'query_name' in annotation_data.columns:
            # Check if the gene IDs are numeric, if yes add 'gene_' in front of them.
            numeric_row_dataframe = pd.to_numeric(annotation_data['query_name'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query_name'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query_name']
            annotation_dict = annotation_data.set_index('query_name')[to_extract_annotations].to_dict('index')
        # 'query' added for compatibility with eggnog-mapper 2.1.2
        elif 'query' in annotation_data.columns:
            numeric_row_dataframe = pd.to_numeric(annotation_data['query'], errors='coerce').notnull()
            if bool(numeric_row_dataframe.any()) is True:
                annotation_data.loc[numeric_row_dataframe, 'query'] = 'gene_' + annotation_data.loc[numeric_row_dataframe, 'query']
            annotation_dict = annotation_data.set_index('query')[to_extract_annotations].to_dict('index')
        for key in annotation_dict:
            yield key, annotation_dict[key]


def create_eggnog_network(eggnog_folder, output_folder, reference_folder=False):
    """From a folder containing eggnog-mapper annotation files, reconstruct draft metabolic networks.

    Args:
        eggnog_folder (str): path to eggnog-mapper annotation files folder
        output_folder (str): path to the output folder
        reference_folder (str): path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|eggnog| Begin KEGG metabolism mapping for folder {0}.'.format(eggnog_folder))

    if reference_folder is not False:
        kegg_model_path = reference_folder
        logger.info('|kegg2bipartitegraph|eggnog| Use reference KEGG model given at {0}.'.format(reference_folder))
    else:
        logger.info('|kegg2bipartitegraph|eggnog| Use default reference KEGG model from {0}.'.format(DATA_ROOT))
        kegg_model_path = os.path.join(DATA_ROOT, 'kegg_model')

    # Download Uniprot metadata and create a json file containing them.
    options = {'eggnog_folder': eggnog_folder, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__

    kegg2bipartitegraph_eggnog_metadata = {}
    kegg2bipartitegraph_eggnog_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    # Check if KEGG model files exist if not create them.
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_eggnog_metadata['reference_path'] = kegg_model_path

    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_json_model_path = os.path.join(kegg_model_path, 'kegg_metadata.json')

    with open(kegg_json_model_path, 'r') as input_metadata_json:
        json_data = json.load(input_metadata_json)

    kegg2bipartitegraph_eggnog_metadata['json_reference_metadata'] = json_data

    kegg_pathways = {}
    with open(kegg_pathways_path, 'r') as open_kegg_pathways_path:
        csvreader = csv.reader(open_kegg_pathways_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            pathway_id = line[0]
            pathway_name = line[1]
            pathway_reactions = line[2].split(',')
            kegg_pathways[pathway_id] = (pathway_name, pathway_reactions)

    kegg_modules = {}
    with open(kegg_modules_path, 'r') as open_kegg_modules_path:
        csvreader = csv.reader(open_kegg_modules_path, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            moudle_id = line[0]
            module_name = line[1]
            module_reactions = line[3].split(',')
            kegg_modules[moudle_id] = (module_name, module_reactions)

    # Create SBML output folder.
    sbml_output_folder_path = os.path.join(output_folder, 'sbml')
    is_valid_dir(sbml_output_folder_path)

    # Create GRAPHML output folder.
    graphml_output_folder_path = os.path.join(output_folder, 'graphml')
    is_valid_dir(graphml_output_folder_path)

    # Create KO output folder.
    ko_output_folder_path = os.path.join(output_folder, 'ko')
    is_valid_dir(ko_output_folder_path)

    # Create pathways annotated output folder.
    pathways_output_folder_path = os.path.join(output_folder, 'pathways')
    is_valid_dir(pathways_output_folder_path)

    # Create modules annotated output folder.
    modules_output_folder_path = os.path.join(output_folder, 'modules')
    is_valid_dir(modules_output_folder_path)

    ko_to_reactions, ec_to_reactions = retrieve_mapping_dictonaries(kegg_rxn_mapping_path)

    # Retrieve EC, KO and add reactions.
    for tsv_file in os.listdir(eggnog_folder):
        enzymes = {}
        org_kos = {}
        taxon_reactions = {}
        if 'emapper.annotations' in tsv_file:
            base_name = tsv_file.replace('.emapper.annotations', '')
            eggnog_annotation_path = os.path.join(eggnog_folder, tsv_file)
            eggnog_annotations = read_annotation(eggnog_annotation_path)
            for gene, annot_dict in eggnog_annotations:
                if 'EC' in annot_dict and annot_dict['EC'] != '-':
                    enzymes[gene] = annot_dict['EC'].split(',')
                if 'KEGG_ko' in annot_dict and annot_dict['KEGG_ko'] != '-':
                    org_kos[gene] = annot_dict['KEGG_ko'].split(',')
                if 'KEGG_Reaction' in annot_dict and annot_dict['KEGG_Reaction'] != '-':
                    for kegg_reaction in annot_dict['KEGG_Reaction'].split(','):
                        if kegg_reaction not in taxon_reactions:
                            taxon_reactions[kegg_reaction] = [gene]
                        else:
                            taxon_reactions[kegg_reaction].append(gene)
            nb_reactions_from_eggnog = len(taxon_reactions)
            logger.info('|kegg2bipartitegraph|eggnog| Added {0} reactions for organism {1} from KEGG reaction of eggnog output.'.format(nb_reactions_from_eggnog, base_name))

            # Use EC found to be associated to reference protein to retrieve KEGG reaction.
            ec_added_reactions = []
            all_ecs = []
            for gene_id in enzymes:
                ec_ids = enzymes[gene_id]
                all_ecs.extend(ec_ids)

                for ec_id in ec_ids:
                    if ec_id in ec_to_reactions:
                        reaction_ids = ec_to_reactions[ec_id]
                        for reaction_id in reaction_ids:
                            ec_added_reactions.append(reaction_id)
                            if reaction_id not in taxon_reactions:
                                taxon_reactions[reaction_id] = [gene_id]
                            else:
                                taxon_reactions[reaction_id].append(gene_id)
            logger.info('|kegg2bipartitegraph|eggnog| Added {0} reactions from {1} EC for organism {2} from eggnog output.'.format(len(set(ec_added_reactions)), len(set(all_ecs)) , base_name))

            # Use KO found in eggnog file to retrieve KEGG reaction.
            ko_added_reactions = []
            all_kos = []
            for gene_id in org_kos:
                ko_ids = org_kos[gene_id]
                all_kos.extend(ko_ids)

                for ko_id in ko_ids:
                    ko_id = ko_id.replace('ko:', '')
                    if ko_id in ko_to_reactions:
                        reaction_ids = ko_to_reactions[ko_id]
                        for reaction_id in reaction_ids:
                            ko_added_reactions.append(reaction_id)
                            if reaction_id not in taxon_reactions:
                                taxon_reactions[reaction_id] = [gene_id]
                            else:
                                taxon_reactions[reaction_id].append(gene_id)
            logger.info('|kegg2bipartitegraph|eggnog| Added {0} reactions from {1} KO for organism {2} from eggnog output.'.format(len(set(ko_added_reactions)), len(set(all_kos)) , base_name))

            total_added_reactions = list(taxon_reactions.keys())
            logger.info('|kegg2bipartitegraph|eggnog| A total of {0} unique reactions are added from EC and KO for taxon {1}.'.format(len(total_added_reactions), base_name))

            # Create pathway file contening pathway with reactions in the taxon.
            pathways_output_file_path = os.path.join(pathways_output_folder_path, base_name+'.tsv')
            organism_pathways = write_pathway_file(kegg_pathways, pathways_output_file_path, total_added_reactions)

            # Create module file contening module with reactions in the taxon.
            modules_output_file_path = os.path.join(modules_output_folder_path, base_name+'.tsv')
            organism_modules = write_module_file(kegg_modules, modules_output_file_path, total_added_reactions)

            kegg_document, kegg_model = create_sbml_from_kegg_reactions(kegg_sbml_model_path, taxon_reactions, organism_pathways, organism_modules)

            # Create file if there is at least 1 reaction.
            if len(kegg_model.getListOfReactions()) > 0:
                # Create SBML file.
                sbml_output_file_path = os.path.join(sbml_output_folder_path, base_name+'.sbml')
                graphml_output_file_path = os.path.join(graphml_output_folder_path, base_name+'.graphml')
                libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
                sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
                logger.info('|kegg2bipartitegraph|eggnog| Network of {0} contains {1} reactions and {2} metabolites.'.format(base_name, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
            else:
                logger.info('|kegg2bipartitegraph|eggnog| No reactions in model for {0}, no SBML file will be created.'.format(base_name))

    kegg_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    compute_stat_kegg(sbml_output_folder_path, kegg_stat_file)

    endtime = time.time()

    duration = endtime - starttime

    kegg2bipartitegraph_eggnog_metadata['kegg2bipartitegraph_eggnog_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_eggnog_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_eggnog_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|eggnog| Draft networks creation complete.')
