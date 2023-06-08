# Copyright (C) 2021-2023 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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

from kegg2bipartitegraph.reconstruct_from_esmecata import create_draft_networks
from kegg2bipartitegraph.reference import create_reference_base
from kegg2bipartitegraph.utils import is_valid_dir
from kegg2bipartitegraph import __version__ as VERSION

MESSAGE = '''
Reconstruct draft metabolic networks using KEGG.
'''
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def main():
    start_time = time.time()

    parser = argparse.ArgumentParser(
        'kegg2bipartitegraph',
        description=MESSAGE + ' For specific help on each subcommand use: esmecata {cmd} --help'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s ' + VERSION + '\n')


    parent_parser_i = argparse.ArgumentParser(add_help=False)
    parent_parser_i.add_argument(
        '-i',
        '--input',
        dest='input',
        required=True,
        help='Input folder corresponds to esmecata annotation folder.',
        metavar='INPUT_DIR')

    parent_parser_o = argparse.ArgumentParser(add_help=False)
    parent_parser_o.add_argument(
        '-o',
        '--output',
        dest='output',
        required=True,
        help='Output directory path.',
        metavar='OUPUT_DIR')

    parent_parser_map_ko = argparse.ArgumentParser(add_help=False)
    parent_parser_map_ko.add_argument(
        '--map-ko',
        dest='map_ko',
        help='From UniProt protein ID, retrieve KEGG ortholgos to infer KEGG reaction from them.',
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

    # subparsers
    subparsers = parser.add_subparsers(
        title='subcommands',
        description='valid subcommands:',
        dest='cmd')

    kegg_reference_parser = subparsers.add_parser(
        'reference',
        help='Create reference data from KEGG.',
        parents=[],
        allow_abbrev=False)

    kegg_esmecata_parser = subparsers.add_parser(
        'esmecata',
        help='Create networks from esmecata results.',
        parents=[
            parent_parser_i, parent_parser_o, parent_parser_map_ko,
            parent_parser_recreate_kegg
            ],
        allow_abbrev=False)

    args = parser.parse_args()

    # If no argument print the help.
    if len(sys.argv) == 1 or len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)

    if args.cmd in ['esmecata']:
        is_valid_dir(args.output)

    formatter = logging.Formatter('%(message)s')
    if args.cmd in ['esmecata']:
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

    if args.cmd == 'esmecata':
        create_draft_networks(args.input, args.output, args.map_ko, args.recreate_kegg)

    if args.cmd == 'reference':
        create_reference_base()

    logger.info("--- Total runtime %.2f seconds ---" % (time.time() - start_time))
    if args.cmd in ['esmecata']:
        logger.warning(f'--- Logs written in {log_file_path} ---')


if __name__ == '__main__':
    main()
