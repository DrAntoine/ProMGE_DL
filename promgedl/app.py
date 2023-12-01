#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import logging


from promgedl.version import __version__
from promgedl.download import download
from promgedl.misc import configureHandler



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

    parser_check = subparser.add_parser(
        "check", 
        prog = f"{softwareName} check",
        description= "local database check module (Sequences)",
        help=f"{softwareName}: check module"
    )

    #region downloadParser
    parser_download.add_argument(
        "-o",
        "--output",
        type=str,
        help="database save location",
        required=True
    )

    parser_download.add_argument(
        "-r",
        "--resume",
        action="store_true",
        help="Allows you to resume downloading. If it is not used, the program downloads all the sequences, even those already present on the disk"
    )

    parser_download.add_argument(
        "--apikey",
        type=str,
        default="",
        help="The NCBI API key allows you to go from 3 requests per second to 10. Please note, the email and API key pair must be linked to your NCBI account otherwise the server may block you. It is not useful to put this parameter if you do not have an API key, the server will block you at 3 requests per second"
    )

    parser_download.add_argument(
        "-@",
        "--email",
        type=str,
        required=True,
        help="email adress to use Entrez API"
    )

    parser_download.add_argument(
        "-p",
        "--promgeURL",
        type=str,
        default="https://promge.embl.de/data/",
        help="URL of the promge site where the files to download are located",
        required=False
    )

    parser_download.add_argument(
        "--distantFiles",
        nargs="+",
        default=["mges.txt.zip", "nested_is-tn_integrons.txt.zip", "arg_mge.txt.zip", "HGT_data.txt.zip", "taxonomy.txt.zip"],
        help="Lists of files to download to the address indicated by the -p argument"
    )

    parser_download.add_argument(
        "--onlyIndex",
        action="store_true",
        default=False,
        help="Prepares the database destination folder and downloads the index file if necessary. The program then stops",
        required=False
    )

    parser_download.add_argument(
        "--showSeqID",
        action="store_true",
        default=False,
        help="Displays the sequence identifiers that will be downloaded and stops the program",
        required=False
    )

    parser_download.add_argument(
        "--seqID",
        type=str,
        default=["*"],
        help="Only sequences with the identifiers specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--genomeID",
        type=str,
        default=["*"],
        help="Only sequences with the genome identifiers specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--taxID",
        type=str,
        default=["*"],
        help="Only sequences with the taxonomy identifiers specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--mgeMinLen",
        type=int,
        help="Only sequences with at least one mge whose size is strictly greater than the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--mgeMaxLen",
        type=int,
        help="Only sequences with at least one mge whose size is less than or equal to the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--minNbProt",
        type=int,
        help="Only sequences containing a number of proteins strictly greater than the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--maxNbProt",
        type=int,
        help="Only sequences with a number of proteins less than or equal to the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--minNbPhageGene",
        type=int,
        help="Only sequences containing a number of phage genes strictly greater than the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--maxNbPhageGene",
        type=int,
        help="Only sequences with a number of phage genes less than or equal to the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--minNbConjGene",
        type=int,
        help="Only sequences containing a number of conjugation genes strictly greater than the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--maxNbConjGene",
        type=int,
        help="Only sequences with a number of conjugation genes less than or equal to the defined value will be downloaded."
    )
    parser_download.add_argument(
        "--minCount",
        type=int,
        help="Description in progress"
    )
    parser_download.add_argument(
        "--maxCount",
        type=int,
        help="Description in progress"
    )
    parser_download.add_argument(
        "--recombinase",
        type=str,
        default=["*"],
        help="Only sequences with the recombinase identifiers specified here will be downloaded if they are present in the index file\nThis parameter is based on an exact match of recombinase identifiers present in the index file. We recommend first downloading the index with the --onlyIndex parameter and then specifying the recombinase identifiers",
        nargs="+"
    )
    parser_download.add_argument(
        "--mgeCategory",
        type=str,
        default=["*"],
        help="Only sequences with the mge categories specified here will be downloaded if they are present in the index file\nThis parameter is based on an exact match of mge categories present in the index file. We recommend first downloading the index with the --onlyIndex parameter and then specifying the mge categories",
        nargs="+"
    )
    parser_download.add_argument(
        "--specI",
        type=str,
        default=["*"],
        help="Only sequences with the specI identifiers specified here will be downloaded if they are present in the index file\nThis parameter is based on an exact match of mge categories present in the index file. We recommend first downloading the index with the --onlyIndex parameter and then specifying the specI identifiers",
        nargs="+"
    )
    parser_download.add_argument(
        "--kingdomLevel",
        type=str,
        default=["*"],
        help="Only sequences in the kingdom(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--phylumLevel",
        type=str,
        default=["*"],
        help="Only sequences in the phylum(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--classLevel",
        type=str,
        default=["*"],
        help="Only sequences in the class(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--orderLevel",
        type=str,
        default=["*"],
        help="Only sequences in the order(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--familyLevel",
        type=str,
        default=["*"],
        help="Only sequences in the family(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--genusLevel",
        type=str,
        default=["*"],
        help="Only sequences in the genus(s) specified here will be downloaded if they are present in the index file",
        nargs="+"
    )
    parser_download.add_argument(
        "--speciesLevel",
        type=str,
        default=["*"],
        help="Only sequences in the species(s) specified here will be downloaded if they are present in the index file.",
        nargs="+"
    )
    #endregion


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
            logging.basicConfig(level=logging.ERROR, handlers=[configureHandler()])
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG, handlers=[configureHandler()])
        else:
            logging.basicConfig(level=logging.INFO, handlers=[configureHandler()])

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