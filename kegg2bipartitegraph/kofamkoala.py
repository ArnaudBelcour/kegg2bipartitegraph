# Copyright (C) 2021-2023 Arnaud Belcour - Inria Dyliss
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
import pandas as pd
import sys
import libsbml

from kegg2bipartitegraph.utils import is_valid_dir
from kegg2bipartitegraph import __version__ as kegg2bipartitegraph_version
from kegg2bipartitegraph.mapping import retrieve_mapping_dictonaries, compute_stat_kegg
from kegg2bipartitegraph.reference import sbml_to_graphml
from kegg2bipartitegraph.sbml import create_sbml_from_kegg_reactions

logger = logging.getLogger(__name__)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_ARCHIVE = os.path.join(*[ROOT, 'data', 'kegg_model.zip'])


def read_kofam_koala_txt(kofam_koala_result_file):
    """Read a kofam koala annotation file and retrieve KO with gene ID.
    Very simple, I think it could easily break with other result files.

    Args:
        kofam_koala_result_file (str): path to kofam koala annotation file
    Returns:
        dict: dict of genes and their annotations as {gene1:[KO1, KO2]}
    """
    df = pd.read_fwf(kofam_koala_result_file, comment='#', sep=' ', header=None)
    # In my example, columns 1 and 2 correspond to gene_id and ko but should be tested on more result files.
    df = df[[1,2]]
    df.columns = ['gene_id', 'ko']
    gene_to_kos = df.groupby('gene_id')['ko'].apply(list).to_dict()

    return gene_to_kos


def create_kofamkoala_network(kofam_koala_folder, output_folder):
    """From a folder containing kofam koala result files, reconstruct draft metabolic networks.

    Args:
        kofam_koala_folder (str): path to kofam koala annotation files folder
        output_folder (str): path to the output folder
    """
    starttime = time.time()
    logger.info('|kegg2bipartitegraph|kofamkoala| Begin KEGG metabolism mapping for folder {0}.'.format(kofam_koala_folder))

    # Download Uniprot metadata and create a json file containing them.
    options = {'kofam_koala_folder': kofam_koala_folder, 'output_folder': output_folder}

    options['tool_dependencies'] = {}
    options['tool_dependencies']['python_package'] = {}
    options['tool_dependencies']['python_package']['Python_version'] = sys.version
    options['tool_dependencies']['python_package']['kegg2bipartitegraph'] = kegg2bipartitegraph_version
    options['tool_dependencies']['python_package']['libsbml'] = libsbml.__version__

    kegg2bipartitegraph_kofamkoala_metadata = {}
    kegg2bipartitegraph_kofamkoala_metadata['tool_options'] = options
    is_valid_dir(output_folder)

    data_kegg_model_path = DATA_ROOT
    # Check if KEGG model files exist if not create them.
    kegg_model_path = os.path.join(data_kegg_model_path, 'kegg_model')
    is_valid_dir(kegg_model_path)
    kegg2bipartitegraph_kofamkoala_metadata['reference_path'] = data_kegg_model_path

    compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    kegg_sbml_model_path = os.path.join(kegg_model_path, 'kegg_model.sbml')
    kegg_rxn_mapping_path = os.path.join(kegg_model_path, 'kegg_mapping.tsv')
    kegg_pathways_path = os.path.join(kegg_model_path, 'kegg_pathways.tsv')
    kegg_modules_path = os.path.join(kegg_model_path, 'kegg_modules.tsv')
    kegg_json_model_path = os.path.join(kegg_model_path, 'kegg_metadata.json')

    with open(kegg_json_model_path, 'r') as input_metadata_json:
        json_data = json.load(input_metadata_json)

    kegg2bipartitegraph_kofamkoala_metadata['json_reference_metadata'] = json_data

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
    for kofamkoal_result_file in os.listdir(kofam_koala_folder):
        taxon_reactions = {}

        base_name = os.path.splitext(kofamkoal_result_file)[0]
        eggnog_annotation_path = os.path.join(kofam_koala_folder, kofamkoal_result_file)
        gene_to_kos = read_kofam_koala_txt(eggnog_annotation_path)

        # Use KO found to be associated to reference protein to retrieve KEGG reaction.
        ko_added_reactions = []
        all_kos = []
        for gene_id in gene_to_kos:
            ko_ids = gene_to_kos[gene_id]
            all_kos.extend(ko_ids)

            for ko_id in ko_ids:
                if ko_id in ko_to_reactions:
                    reaction_ids = ko_to_reactions[ko_id]
                    for reaction_id in reaction_ids:
                        ko_added_reactions.append(reaction_id)
                        if reaction_id not in taxon_reactions:
                            taxon_reactions[reaction_id] = [gene_id]
                        else:
                            taxon_reactions[reaction_id].append(gene_id)
        logger.info('|kegg2bipartitegraph|kofamkoala| Added {0} reactions from {1} KO for organism {2} from eggnog output.'.format(len(set(ko_added_reactions)), len(set(all_kos)) , base_name))

        total_added_reactions = list(taxon_reactions.keys())

        # Create pathway file contening pathway with reactions in the taxon.
        pathways_output_file_path = os.path.join(pathways_output_folder_path, base_name+'.tsv')
        with open(pathways_output_file_path, 'w') as open_pathways_output_file_path:
            csvwriter = csv.writer(open_pathways_output_file_path, delimiter='\t')
            csvwriter.writerow(['pathway_id', 'pathway_name', 'pathway_completion_ratio', 'pathway_reaction_in_taxon', 'pathway_reaction'])
            for pathway in kegg_pathways:
                pathway_reactions = kegg_pathways[pathway][1]
                pathway_reaction_in_taxon = set(pathway_reactions).intersection(set(total_added_reactions))
                if len(pathway_reaction_in_taxon) > 0:
                    pathway_name = kegg_pathways[pathway][0]
                    pathway_completion_ratio = len(pathway_reaction_in_taxon) / len(pathway_reactions)
                    csvwriter.writerow([pathway, pathway_name, pathway_completion_ratio, ','.join(pathway_reaction_in_taxon), ','.join(pathway_reactions)])

        # Create module file contening module with reactions in the taxon.
        modules_output_file_path = os.path.join(modules_output_folder_path, base_name+'.tsv')
        with open(modules_output_file_path, 'w') as open_modules_output_file_path:
            csvwriter = csv.writer(open_modules_output_file_path, delimiter='\t')
            csvwriter.writerow(['module_id', 'module_name', 'module_completion_ratio', 'module_reaction_in_taxon', 'module_reaction'])
            for module in kegg_modules:
                module_reactions = kegg_modules[module][1]
                module_reaction_in_taxon = set(module_reactions).intersection(set(total_added_reactions))
                if len(module_reaction_in_taxon) > 0:
                    module_name = kegg_modules[module][0]
                    module_completion_ratio = len(module_reaction_in_taxon) / len(module_reactions)
                    csvwriter.writerow([module, module_name, module_completion_ratio, ','.join(module_reaction_in_taxon), ','.join(module_reactions)])

        kegg_document, kegg_model = create_sbml_from_kegg_reactions(kegg_sbml_model_path, taxon_reactions)

        # Create file if there is at least 1 reaction.
        if len(kegg_model.getListOfReactions()) > 0:
            # Create SBML file.
            sbml_output_file_path = os.path.join(sbml_output_folder_path, base_name+'.sbml')
            graphml_output_file_path = os.path.join(graphml_output_folder_path, base_name+'.graphml')
            libsbml.writeSBMLToFile(kegg_document, sbml_output_file_path)
            sbml_to_graphml(sbml_output_file_path, graphml_output_file_path)
            logger.info('|kegg2bipartitegraph|kofamkoala| Network of {0} contains {1} reactions and {2} metabolites.'.format(base_name, len(kegg_model.getListOfReactions()), len(kegg_model.getListOfSpecies())))
        else:
            logger.info('|kegg2bipartitegraph|kofamkoala| No reactions in model for {0}, no SBML file will be created.'.format(base_name))


    clust_stat_file = os.path.join(output_folder, 'stat_number_kegg.tsv')
    compute_stat_kegg(sbml_output_folder_path, clust_stat_file)

    endtime = time.time()

    duration = endtime - starttime

    kegg2bipartitegraph_kofamkoala_metadata['kegg2bipartitegraph_kofamkoala_duration'] = duration
    uniprot_metadata_file = os.path.join(output_folder, 'kegg2bipartitegraph_kofamkoala_kegg.json')
    with open(uniprot_metadata_file, 'w') as ouput_file:
        json.dump(kegg2bipartitegraph_kofamkoala_metadata, ouput_file, indent=4)
    logger.info('|kegg2bipartitegraph|kofamkoala| Draft networks creation complete.')
