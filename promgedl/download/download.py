import os, logging, signal

from promgedl.download.tools import * 
logger = logging.getLogger(__name__)


def handler(signum, frame):
    logger.critical("The stop signal has just been used!")
    exit(1)

signal.signal(signal.SIGINT, handler)


def run(args):
    """
    main function for promgedl download
    """
    
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
        logger.info("Creating database output folder")
        for ds in directory_structure:
            logger.debug(f"Creating {directory_structure[ds]}")
            os.makedirs(directory_structure[ds])

    for file in args.distantFiles:
        if not checkIndexFile(f"{directory_structure['indexFolder']}/{file}"):
            downloadIndex(baseURL = args.promgeURL, distantFile=file, indexLocation=f"{directory_structure['indexFolder']}")
    
    if args.onlyIndex:
        exit()
    sequencesToDownload = []
    if 'mges.txt.zip' in args.distantFiles:
        sequencesToDownload = extractSequencesIdFromIndexFile(f"{directory_structure['indexFolder']}/mges.txt", args)

    if len(sequencesToDownload) == 0:
        logger.error("The filters applied are too strict! No sequence was accepted")
        exit()
    sequencesToDownloadFiltered = extractSequenceIDFromTree(makeFilteredTree(sequencesToDownload))
    sequencesToDownloadFiltered.sort()
    # Extraire les n° de seq (make/extract filtered tree)
    if args.exportSeqID:
        path = f"directory_structure['indexFolder']/{datetime.date.today()}_seqIDsExported"
        if args.exportSeqID != "":
            path = args.exportSeqID
        exportSeqIDs(sequencesToDownloadFiltered, path, args)
    
    if args.showSeqID:
        displaySeqIds(sequencesToDownloadFiltered)
        exit()

    downloadSeq(sequencesToDownloadFiltered, f"{directory_structure['sequencesFolder']}/", f"{directory_structure['indexFolder']}",args.email, args.resume, entrezAPI=args.apikey)
    logger.info("Goodbye")