#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os, sys
import logging


from promgedl.version import __version__

def main():
    description = "Script to download the ProMGE database"
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

    parser.add_argument(
        "-p",
        "--promge",
        type=str,
        default="https://promge.embl.de/data/mges.txt.zip",
        help="url of mges index file",
        required=False
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="database save location",
        required=True
    )

    parser.add_argument(
        "-@",
        "--email",
        type=str,
        required=True,
        help="email adress to use Entrez API"
    )
    parserlogging = parser.add_mutually_exclusive_group()
    parserlogging.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        default=False,
        help="disable info logging"
    )

    parserlogging.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="overwhelms the terminal with debug information"
    )

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
        logging.shutdown()
    except AttributeError as err:
        logger = logging.getLogger(__name__)
        logger.debug(err)
        parser.print_help()
        logging.shutdown()

if __name__ == "__main__":
    main()