import os, logging

from Bio import Entrez

from promgedl.download.tools import check_directory_structure

def run(args):
    """
    main function for promgedl download
    """
    
    logger = logging.getLogger(__name__)
    logger.debug("Welcome in promgedl download")

    logger.debug("Checking if the output folder already exist")

    directory_structure = [
        f"{args.output}",
        f"{args.output}/index",
        f"{args.output}/sequences"
    ]

    structureHealth = check_directory_structure(directory_structure)
    if structureHealth == "Complete":
        logger.warning(f"{args.output} already exist!")
    elif structureHealth == "Partial":
        logger.error(f"{args.output} exist but is incomplete ! Please remove this folder or change the destination path to avoid unpredictable behaviors")
        exit(1)
    elif structureHealth == "Absent":
        logger.info(f"Creating database output folder")
        for ds in directory_structure:
            logger.debug(f"Creating {ds}")
            os.makedirs(ds)

    # Télécharger le fichier proMGE
    # Unzip le fichier proMGE
    # Extraire les n° de seq (make/extract filtered tree)
    

    
    pass