# Copyright (C) 2021-2025 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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
import logging
import os
import sys
import time

from kegg2bipartitegraph.esmecata import create_esmecata_network
from kegg2bipartitegraph.organism import create_organism_network
from kegg2bipartitegraph.eggnog import create_eggnog_network
from kegg2bipartitegraph.kofamkoala import create_kofamkoala_network
from kegg2bipartitegraph.picrust import create_picrust_network
from kegg2bipartitegraph.genbank import create_gbff_network
from kegg2bipartitegraph.reference import create_reference_base
from kegg2bipartitegraph.scope import compute_scope
from kegg2bipartitegraph.utils import is_valid_dir, get_internal_reference_model_path, get_current_database_version
from kegg2bipartitegraph import __version__ as VERSION

MESSAGE = '''
Reconstruct draft metabolic networks using KEGG.
'''
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

ROOT = os.path.dirname(__file__)
DATA_ROOT = os.path.join(ROOT, 'data')
KEGG_MODEL = os.path.join(DATA_ROOT, 'kegg_model')
SEED_FILES = {'HYPERTHERM_MED': os.path.join(KEGG_MODEL, 'seed_hyperterm_med.txt'), 'MESO_MEDIUM': os.path.join(KEGG_MODEL, 'seeds_meso_medium.txt'),
              'PSYCHROMED': os.path.join(KEGG_MODEL, 'seeds_psychromed.txt'), 'THERM_MED': os.path.join(KEGG_MODEL, 'seeds_therm_med.txt'),
              'UNION': os.path.join(KEGG_MODEL, 'seed_union_all.txt')}

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'k2bg',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help'
    )

    # Retrieve version of KEGG used for reference model.
    reference_model_path = get_internal_reference_model_path()
    database_version = get_current_database_version(reference_model_path)
    if database_version is not False:
        kegg_database_version = 'Version of KEGG database used for reference model: ' + database_version
        version_message = '%(prog)s version: ' + VERSION + ', '+ kegg_database_version
    else:
        version_message = '%(prog)s version: ' + VERSION + ', no database version found'

    parser.add_argument(
        '--version',
        action='version',
        version=version_message)

    parent_parser_i = argparse.ArgumentParser(add_help=False)
    parent_parser_i.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input folder corresponding to esmecata annotation folder.',
        metavar='INPUT_DIR')

    parent_parser_i_scope = argparse.ArgumentParser(add_help=False)
    parent_parser_i_scope.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input folder containing graphml file.',
        metavar='INPUT_DIR')

    parent_parser_s_scope = argparse.ArgumentParser(add_help=False)
    parent_parser_s_scope.add_argument(
        '-s',
        '--seed',
        dest='seed_file',
        required=False,
        help='Seed file for scope analysis (txt file of KEGG IDs, one per row). By default, use ones from kegg2bipartitegraph. If you want a specific ones, select one among: MESO_MEDIUM, PSYCHROMED, HYPERTHERM_MED, THERM_MED',
        metavar='INPUT_FILE')

    parent_parser_i_eggnog = argparse.ArgumentParser(add_help=False)
    parent_parser_i_eggnog.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input folder corresponding to a folder containing eggnog-mapper .emapper.annotations files.',
        metavar='INPUT_DIR')

    parent_parser_i_kofamkoala = argparse.ArgumentParser(add_help=False)
    parent_parser_i_kofamkoala.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input folder corresponding to a folder containing multiple kofam koala result files.',
        metavar='INPUT_DIR')

    parent_parser_i_org = argparse.ArgumentParser(add_help=False)
    parent_parser_i_org.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='KEGG organism code (for example hsa for Homo sapiens (human)).',
        metavar='INPUT_STR')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_o_reference = argparse.ArgumentParser(add_help=False)
    parent_parser_o_reference.add_argument(
        '-o',
        '--output',
        dest='output',
        required=False,
        help='Redirect the creation of the reference database to another folder from the one inside the package repository (default: None).',
        metavar='OUPUT_DIR')

    parent_parser_map_ko = argparse.ArgumentParser(add_help=False)
    parent_parser_map_ko.add_argument(
        '--map-ko',
        dest='map_ko',
        help='From UniProt protein ID, retrieve KEGG ortholgos to infer KEGG reactions.',
        required=False,
        action='store_true',
        default=None)

    parent_parser_recreate_kegg = argparse.ArgumentParser(add_help=False)
    parent_parser_recreate_kegg.add_argument(
        '--recreate-kegg',
        dest='recreate_kegg',
        help='Recreate KEGG model by downloading and extracting informations from KEGG files.',
        required=False,
        action='store_true',
        default=None)

    parent_parser_r = argparse.ArgumentParser(add_help=False)
    parent_parser_r.add_argument(
        '-r',
        '--reference',
        dest='reference_folder',
        required=False,
        help='Path to a reference KEGG folder, to use it instead of the default ones contained in kegg2bipartitegraph.',
        metavar='REFERENCE_PATH',
        default=False)

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    kegg_reference_parser = subparsers.add_parser(
        'reference',
        help='Create reference data from KEGG.',
        parents=[parent_parser_o_reference],
        allow_abbrev=False)

    kegg_esmecata_parser = subparsers.add_parser(
        'reconstruct_from_esmecata',
        help='Create networks from esmecata results.',
        parents=[
            parent_parser_i, parent_parser_o, parent_parser_map_ko,
            parent_parser_r
            ],
        allow_abbrev=False)

    kegg_organism_parser = subparsers.add_parser(
        'reconstruct_from_organism',
        help='Create network from a KEGG organism code.',
        parents=[
            parent_parser_i_org, parent_parser_o, parent_parser_r
            ],
        allow_abbrev=False)

    kegg_eggnog_parser = subparsers.add_parser(
        'reconstruct_from_eggnog',
        help='Create network from a folder containing multiple eggnog-mapper annotation files.',
        parents=[
            parent_parser_i_eggnog, parent_parser_o, parent_parser_r
            ],
        allow_abbrev=False)

    kegg_kofamkoala_parser = subparsers.add_parser(
        'reconstruct_from_kofamkoala',
        help='Create network from a folder containing multiple kofam koala result files.',
        parents=[
            parent_parser_i_kofamkoala, parent_parser_o, parent_parser_r
            ],
        allow_abbrev=False)

    kegg_picrust_parser = subparsers.add_parser(
        'reconstruct_from_picrust',
        help='Create networks from picrust results.',
        parents=[
            parent_parser_i, parent_parser_o, parent_parser_r
            ],
        allow_abbrev=False)

    kegg_genbank_parser = subparsers.add_parser(
        'reconstruct_from_genbank',
        help='Create networks from GenBank files.',
        parents=[
            parent_parser_i, parent_parser_o, parent_parser_r
            ],
        allow_abbrev=False)

    kegg_scope_parser = subparsers.add_parser(
        'scope',
        help='Compute scope from folder of graphml files.',
        parents=[parent_parser_i_scope, parent_parser_o, parent_parser_s_scope],
        allow_abbrev=False)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    if args.cmd in ['reconstruct_from_esmecata', 'reconstruct_from_organism', 'reconstruct_from_eggnogg',
                    'reconstruct_from_kofamkoala', 'reconstruct_from_picrust', 'reconstruct_from_genbank',
                    'scope']:
        is_valid_dir(args.output)

    formatter = logging.Formatter('%(message)s')
    if args.cmd in ['reconstruct_from_esmecata', 'reconstruct_from_organism', 'reconstruct_from_eggnogg',
                    'reconstruct_from_kofamkoala', 'reconstruct_from_picrust', 'reconstruct_from_genbank',
                    'scope']:
        # add logger in file
        log_file_path = os.path.join(args.output, f'kegg2bipartitegraph{args.cmd}.log')
        file_handler = logging.FileHandler(log_file_path, 'w+')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    # set up the default console logger
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if args.cmd == 'reference':
        create_reference_base(args.output)
    elif args.cmd == 'reconstruct_from_esmecata':
        create_esmecata_network(args.input, args.output, args.map_ko, args.reference_folder)
    elif args.cmd == 'reconstruct_from_organism':
        create_organism_network(args.input, args.output, args.reference_folder)
    elif args.cmd == 'reconstruct_from_eggnog':
        create_eggnog_network(args.input, args.output, args.reference_folder)
    elif args.cmd == 'reconstruct_from_kofamkoala':
        create_kofamkoala_network(args.input, args.output, args.reference_folder)
    elif args.cmd == 'reconstruct_from_picrust':
        create_picrust_network(args.input, args.output, args.reference_folder)
    elif args.cmd == 'reconstruct_from_genbank':
        create_gbff_network(args.input, args.output, args.reference_folder)
    elif args.cmd == 'scope':
        output_folder = args.output
        is_valid_dir(output_folder)
        seed_file = args.seed_file
        input_seed_files = []
        if seed_file is None:
            for seed_file in ['MESO_MEDIUM', 'PSYCHROMED', 'HYPERTHERM_MED', 'THERM_MED', 'UNION']:
                input_seed_files.append(SEED_FILES[seed_file])
        elif ',' in seed_file:
            seed_files = seed_file.split(',')
            for seed_file in seed_files:
                if seed_file in ['MESO_MEDIUM', 'PSYCHROMED', 'HYPERTHERM_MED', 'THERM_MED', 'UNION']:
                    input_seed_files.append(SEED_FILES[seed_file])
                else:
                    input_seed_files.append(seed_file)
        else:
            if seed_file in ['MESO_MEDIUM', 'PSYCHROMED', 'HYPERTHERM_MED', 'THERM_MED', 'UNION']:
                input_seed_files.append(SEED_FILES[seed_file])
            else:
                input_seed_files.append(seed_file)
        for seed_file in input_seed_files:
            seed_name = os.path.splitext(os.path.basename(seed_file))[0]
            seed_output_folder = os.path.join(output_folder, seed_name)
            compute_scope(args.input, seed_output_folder, seed_file)

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    if args.cmd in ['esmecata']:
        logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
