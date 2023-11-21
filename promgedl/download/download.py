import os, logging

from Bio import Entrez

from promgedl.download.tools import * #checkDirectoryStructure, checkIndexFile, downloadIndex, extractSeqencesIdFromIndexFile

def run(args):
    """
    main function for promgedl download
    """
    
    logger = logging.getLogger(__name__)
    logger.debug("Welcome in promgedl download")

    logger.debug("Checking if the output folder already exist")

    directory_structure = {
        "base" : f"{args.output}",
        "indexFolder": f"{args.output}/index",
        "sequencesFolder" :f"{args.output}/sequences"
    }

    structureHealth = checkDirectoryStructure(directory_structure)
    logger.debug(f"StructureHealth : {structureHealth}")

    if structureHealth == "Complete":
        logger.warning(f"{args.output} already exist!")
    elif structureHealth == "Partial":
        logger.error(f"{args.output} exist but is incomplete ! Please remove this folder or change the destination path to avoid unpredictable behaviors")
        exit(1)
    elif structureHealth == "Absent":
        logger.info(f"Creating database output folder")
        for ds in directory_structure.keys():
            logger.debug(f"Creating {directory_structure[ds]}")
            os.makedirs(directory_structure[ds])

    if not checkIndexFile(f"{directory_structure['indexFolder']}/mges.zip"):
        downloadIndex(url = args.promgeURL, indexLocation=f"{directory_structure['indexFolder']}/")
    
    if args.onlyIndex:
        exit()

    sequencesToDownload = extractSequencesIdFromIndexFile(f"{directory_structure['indexFolder']}/mges.txt", args)
    if len(sequencesToDownload) == 0:
        logger.error("The filters applied are too strict! No sequence was accepted")
        exit()
    sequencesToDownloadFiltered = extractSequenceIDFromTree(makeFilteredTree(sequencesToDownload))
    sequencesToDownloadFiltered.sort()
    # Extraire les n° de seq (make/extract filtered tree)
    if args.showSeqID:
        displaySeqIds(sequencesToDownloadFiltered)
        exit()

    downloadSeq(sequencesToDownloadFiltered, f"{directory_structure['sequencesFolder']}/", args.email)
    logger.info("Goodbye")