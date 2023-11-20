#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os, sys
import logging


from promgedl.version import __version__
from promgedl.download import download

def main():
    description = "Script faciliting the use of the ProMGE database"
    softwareName = "promgedl"
    parser = argparse.ArgumentParser(
        prog=softwareName,
        description=description,
        usage=f"{softwareName} [option] ...",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
        required=False
    )
    subparser = parser.add_subparsers(
        title="available subcommands",
        metavar=""
    )

    parser_download = subparser.add_parser(
        "download",
        prog = f"{softwareName} download",
        description = "index and sequence download module",
        help=f"{softwareName}: download module"
    )

    parser_access = subparser.add_parser(
        "access", 
        prog = f"{softwareName} access",
        description= "local database access module (Sequences + MGEs)",
        help=f"{softwareName}: access module"
    )

    parser_download.add_argument(
        "-p",
        "--promge",
        type=str,
        default="https://promge.embl.de/data/mges.txt.zip",
        help="url of mges index file",
        required=False
    )

    parser_download.add_argument(
        "-o",
        "--output",
        type=str,
        help="database save location",
        required=True
    )

    parser_download.add_argument(
        "-@",
        "--email",
        type=str,
        required=True,
        help="email adress to use Entrez API"
    )

    parser_download.add_argument(
        "-u",
        "--update",
        action="store_true",
        default=False,
        help="\033[91m[DANGER]\033[0m Overwrites existing data in the destination folder"
    )

    parser_download.add_argument(
        "-c",
        "--clean",
        action="store_true",
        default=False,
        help="\033[91m[DANGER]\033[0m Overwrites existing data in the destination folder"
    )

    parserloggingDL = parser_download.add_mutually_exclusive_group()
    parserloggingDL.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="disable info logging"
    )

    parserloggingDL.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="overwhelms the terminal with debug information"
    )

    parser_download.set_defaults(func=download.run)

    args = parser.parse_args()

    try:
        if args.version:
            print(f"{__name__} version {__version__}")
            sys.exit(0)
        elif args.quiet:
            logging.basicConfig(level=logging.ERROR)
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        logger.info(f"Using {__name__} version {__version__}")
        logger.debug(f"Using Debug logger")
        args.func(args)
        logging.shutdown()
    except AttributeError as err:
        logger = logging.getLogger(__name__)
        logger.debug(err)
        parser.print_help()
        logging.shutdown()

if __name__ == "__main__":
    main()