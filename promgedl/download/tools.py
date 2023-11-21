import os, requests, shutil, datetime, logging
from zipfile import ZipFile

from tqdm.auto import tqdm
from Bio import Entrez


logger = logging.getLogger(__name__)

class noeud:

    def __init__(self, IdPere, value, idFilsG=None, idFrereD=None) -> None:
        self.pere = IdPere
        self.content=value
        self.filsG = idFilsG
        self.frereD = idFrereD

class DownloadError(Exception):
    pass

def checkDirectoryStructure(structure):
    logger.debug("checkDirectoryStructure function")
    folderStatus = []
    for s in structure.keys():
        folderStatus.append(os.path.isdir(structure[s]))
    if len(set(folderStatus))==2:
        return "Partial"
    else:
        if list(set(folderStatus))[0]==True:
            return "Complete"
        else:
            return "Absent"
    
def downloadIndex(url, indexLocation):
    logger.debug("downloadIndex function")
    with requests.get(url, stream=True) as req:
        totalLength = int(req.headers.get("Content-Length"))

        with tqdm.wrapattr(req.raw, "read", total=totalLength, desc="Downloading index file") as rawdata:
            with open(f"{indexLocation}/mges.zip", "wb") as outputfile:
                shutil.copyfileobj(rawdata, outputfile)
    with ZipFile(f"{indexLocation}/mges.zip","r") as zipObj:
        zipObj.extractall(indexLocation)

def checkIndexFile(indexLocation):
    logger.debug("checkIndexFile function")
    numberOfDaysBeforeDownloadAgain = 7
    logger.debug(indexLocation)
    if not os.path.exists(indexLocation):
        logger.debug(f"The archive does not exit. It will be downloaded")
        return False
    logger.debug("mges zipfile exist")
    lastModified = os.path.getmtime(indexLocation)
    numberofDaysElapsed = datetime.datetime.now() - datetime.datetime.fromtimestamp(lastModified)
    days = numberofDaysElapsed.days
    logger.debug(f"Zipfile downloaded {days} day(s) ago")
    if days >= numberOfDaysBeforeDownloadAgain:
        logger.debug(f"The archive is to old and will be downloaded again. (file's days >= {numberOfDaysBeforeDownloadAgain})")
        return False
    logger.debug(f"The archive is sufficiently recent")
    return True
    
def __count_generator__(reader):
    b = reader(1024*1024)
    while b:
        yield b
        b = reader(1024*1024)  

def estimateFileLength(filepath):
    logger.debug("estimateFileLength function")
    with open(filepath, "rb") as file:
        countGenerator = __count_generator__(file.raw.read)
        linecount = sum(buffer.count(b'\n') for buffer in countGenerator)
    return linecount

def convertToInt(string):
    try:
        return int(string)
    except:
        return -1

def applyFilters(line, dbindex, args):
    mgeSequenceIdentifier = line[0]
    sequenceIdentifier, mgePositionStartEnd = mgeSequenceIdentifier.split(":")
    taxid, genome, seqID = sequenceIdentifier.split(".")
    mge_length = convertToInt(line[dbindex.index("mge_length")])
    number_of_proteins = convertToInt(line[dbindex.index("number_of_proteins")])
    nb_phage_genes = convertToInt(line[dbindex.index("phage_genes")])
    nb_conjugation_genes = convertToInt(line[dbindex.index("conjugation_genes")])
    recombinase = line[dbindex.index("recombinase")]
    mge_category = line[dbindex.index("mge_category")]
    count = convertToInt(line[dbindex.index("count")])
    specI = line[dbindex.index("specI")]
    s_kingdom = line[dbindex.index("kingdom")]
    s_phylum = line[dbindex.index("phylum")]
    s_class = line[dbindex.index("class")]
    s_order = line[dbindex.index("order")]
    s_family = line[dbindex.index("family")]
    s_genus = line[dbindex.index("genus")]
    s_species = line[dbindex.index("species")]

    if args.seqID!=["*"] and seqID not in args.seqID:
        logger.debug("seqID filter")
        return None
    if args.genomeID!=["*"] and genome not in args.genomeID:
        logger.debug("genomeID filter")
        return None
    if args.taxID!=["*"] and taxid not in args.taxID:
        logger.debug(f"taxID filter: {taxid} not in {args.taxID}")
        return None
    
    if args.mgeMinLen and mge_length>=0 and args.mgeMinLen >=0 and mge_length <= args.mgeMinLen:
        logger.debug("mgeMinLen filter")
        return None
    if args.mgeMaxLen and mge_length>=0 and args.mgeMaxLen >=0 and mge_length > args.mgeMaxLen:
        logger.debug("mgeMaxLen filter")
        return None
    if args.minNbProt and number_of_proteins>=0 and args.minNbProt >=0 and number_of_proteins <= args.minNbProt:
        logger.debug("minNbProt filter")
        return None
    if args.maxNbProt and number_of_proteins>=0 and args.maxNbProt >=0 and number_of_proteins > args.maxNbProt:
        logger.debug("maxNbProt filter")
        return None
    if args.minNbPhageGene and nb_phage_genes>=0 and args.minNbPhageGene >=0 and nb_phage_genes <= args.minNbPhageGene:
        logger.debug("minNbPhageGene filter")
        return None
    if args.maxNbPhageGene and nb_phage_genes>=0 and args.maxNbPhageGene >=0 and nb_phage_genes > args.maxNbPhageGene:
        logger.debug("maxNbPhageGene filter")
        return None
    if args.minNbConjGene and nb_conjugation_genes>=0 and args.minNbConjGene >=0 and nb_conjugation_genes <= args.minNbConjGene:
        logger.debug("minNbConjGene filter")
        return None
    if args.maxNbConjGene and nb_conjugation_genes>=0 and args.maxNbConjGene >=0 and nb_conjugation_genes > args.maxNbConjGene:
        logger.debug("maxNbConjGene filter")
        return None
    if args.minCount and count>=0 and args.minCount >=0 and count <= args.minCount:
        logger.debug("minCount filter")
        return None
    if args.maxCount and count>=0 and args.maxCount >=0 and count > args.maxCount:
        logger.debug("maxCount filter")
        return None

    if args.recombinase!=["*"]:
        accept = False
        recombinases = recombinase.split(",")
        for r in recombinases:
            if r in args.recombinase:
                accept = True
        if not accept:
            logger.debug("recombinase filter")
            return None
    if args.mgeCategory!=["*"]:
        accept = False
        mge_categories = mge_category.split(",")
        for m in mge_categories:
            if m in args.mgeCategory:
                accept = True
        if not accept:
            logger.debug("mgeCategory filter")
            return None
    if args.specI!=["*"]:
        accept = False
        specIes = specI.split(",")
        for s in specIes:
            if s in args.specI:
                accept = True
        if not accept:
            logger.debug("specI filter")
            return None
    
    if args.kingdomLevel!=["*"] and s_kingdom not in args.kingdomLevel:
        logger.debug("kingdomLevel filter")
        return None
    if args.phylumLevel!=["*"] and s_phylum not in args.phylumLevel:
        logger.debug("phylumLevel filter")
        return None
    if args.classLevel!=["*"] and s_class not in args.classLevel:
        logger.debug("classLevel filter")
        return None
    if args.orderLevel!=["*"] and s_order not in args.orderLevel:
        logger.debug("orderLevel filter")
        return None
    if args.familyLevel!=["*"] and s_family not in args.familyLevel:
        logger.debug("familyLevel filter")
        return None
    if args.genusLevel!=["*"] and s_genus not in args.genusLevel:
        logger.debug("genusLevel filter")
        return None
    if args.speciesLevel!=["*"] and s_species not in args.speciesLevel:
        logger.debug("speciesLevel filter")
        return None

    return seqID

def extractSequencesIdFromIndexFile(indexPath, args):
    logger.debug("extractSeqencesIdFromIndexFile function")
    actualLine = 0
    indexFileIndex = []
    sequencesIDList = []
    with tqdm(total=estimateFileLength(indexPath), desc="Extracting accession numbers from index file") as pbar:
        with open(indexPath, "r") as indexFile:
            for line in indexFile:
                actualLine+=1
                splitedLine = line.strip().split("\t")
                if actualLine==1:
                    indexFileIndex = splitedLine
                    logger.debug(f"IndexFileIndex : {indexFileIndex}")
                else:
                    sequenceID = applyFilters(splitedLine, indexFileIndex, args)
                    if sequenceID:
                        sequencesIDList.append(sequenceID)
                pbar.update()
    return sequencesIDList

def displaySeqIds(seqIDs):
    line = ""
    for i in range(len(seqIDs)):
        line += f"{seqIDs[i]},"
        if (i+1)%10 == 0:
            logger.info(line[:-1])
            line = ""
    if line != "":
        logger.info(line[:-1])
    logger.info(f"Total: {len(seqIDs)} sequence(s)")
    

def makeFilteredTree(sequences):
    tree = [noeud(IdPere=None, value="START")]
    nbseq = len(sequences)
    with tqdm(total=nbseq, desc="Sorting accession numbers (1/2)") as pbar:
        for i in range(nbseq):
            seqID = sequences[i]
            nodeIndex=0
            for char in seqID:
                node = tree[nodeIndex]
                if node.filsG is None:
                    tree.append(noeud(IdPere=nodeIndex, value=char))
                    tree[nodeIndex].filsG=len(tree)-1
                    nodeIndex=tree[nodeIndex].filsG
                else:
                    childNodeindex = node.filsG
                    childNode = tree[childNodeindex]
                    while childNode.content != char :
                        if childNode.frereD is None:
                            tree.append(noeud(IdPere=nodeIndex, value=char))
                            childNode.frereD = len(tree)-1
                        else:
                            childNodeindex = childNode.frereD
                        childNode = tree[childNodeindex]
                    nodeIndex=childNodeindex
            pbar.update()
    return tree

def extractSequenceIDFromTree(tree):
    sequencesID = []
    # with alive_bar(len(tree)-1) as bar:
    with tqdm(total=len(tree)-1, desc="Sorting accession numbers (2/2)") as pbar:
        while len(tree)>1:
            seq = ""
            node = tree[-1]
            while node.content != "START":
                seq+=node.content
                node = tree[node.pere]
            seq = seq[::-1]
            sequencesID.append(seq)
            node = tree[-1]
            indexToRemove = len(tree)-1
            while node.content != "START":
                if node.filsG is None and node.frereD is None:
                    if tree[node.pere].filsG == indexToRemove:
                        tree[node.pere].filsG = None
                        newindexToRemove=node.pere
                        del tree[indexToRemove]
                        indexToRemove = newindexToRemove
                    else:
                        child=tree[tree[node.pere].filsG]
                        while child.frereD != indexToRemove:
                            child=tree[child.frereD]
                        child.frereD = None
                        del tree[indexToRemove]
                    pbar.update()
                node = tree[node.pere]
    return sequencesID

def searchEntrez(query, entrezDB = "nuccore"):
    print(f"Searching Entrez IDs db:{entrezDB}")
    subQueryString = list(query)
    queryString = ", ".join(subQueryString)
    handle = Entrez.esearch(db=entrezDB, rettype="fasta", retmode="text", term=queryString, idtype="acc")
    records = Entrez.read(handle)
    handle.close()
    return records["IdList"]

def fetchEntrez(entrezIDs, entrezDB = "nuccore"):
    print(f"Fetching sequences on db:{entrezDB}")
    handle = Entrez.efetch(db=entrezDB, id=entrezIDs, rettype="fasta", retmode="text")
    data = handle.readlines()
    handle.close()
    return data

def ParseSequences(stream):
    print("Parsing data")
    sequences = {}
    seqID = None
    buffer = []
    for line in stream:
        if ">" in line:
            if seqID:
                sequences[seqID] = buffer
                buffer = []
            seqid = line[1:].split(" ")[0]
        buffer.append(line)
    sequences[seqid] = buffer
    return sequences

def getDestinationFolder(sequence,dbLocation):
    indexLvl1 = sequence[:2]
    indexLvl2 = sequence[:4]
    seqPath=f"{dbLocation}/{indexLvl1}/{indexLvl2}"
    for i in range(2):
        if not os.path.exists(seqPath):
            if i == 0: os.makedirs(seqPath)
            else: raise ValueError(seqPath)
    return seqPath

def writeSequence(seqPath, seqName, sequence):
    with open(f"{seqPath}/{seqName}.fasta", "w") as seqFile:
        for line in sequence:
            seqFile.write(line)

def downloadSeq(seqIDs, dbLocation, entrezEmail, entrezAPI=""):
    downloadByStep = 25
    failedDownload = []
    entrezDB = "nuccore"
    Entrez.email = entrezEmail
    if entrezAPI != "":
        Entrez.api_key = entrezAPI

    with tqdm(total=len(seqIDs), desc="Downloading sequences") as pbar:
        for step in range(0, len(seqIDs), downloadByStep):
            query = seqIDs[step: step+downloadByStep]
            fastaDict = None
            for downloadAttempt in range(3):
                try:
                    entrezIDs = searchEntrez(query, entrezDB)
                    if entrezIDs == []:
                        raise DownloadError("No ID found in Entrez database")
                    rawsequences = fetchEntrez(entrezIDs, entrezDB)
                    if len(rawsequences)<=1:
                        raise DownloadError(f"fetchEntrez Error : {rawsequences[0]}")
                    fastaDict = ParseSequences(rawsequences)
                    break
                except DownloadError as err:
                    logger.error(err)
                logger.info(f"Download attempt {downloadAttempt} failed")
            
            if fastaDict is None:
                failedDownload += query
                pbar.update(len(query))
                logger.warning(f"Unable to download these sequences: {', '.join(query)}")
            else:
                fastaSeqKey = [seqID for seqID in fastaDict.keys()]
                fastaSeqKeyShort = [seqID.split(".")[0] for seqID in fastaSeqKey]
                localFailedDownload = []
                for seqID in query:
                    if seqID not in fastaSeqKeyShort:
                        localFailedDownload.append(seqID)
                    else:
                        seqPath = getDestinationFolder(seqID, dbLocation)
                        writeSequence(seqPath, seqID, fastaDict[fastaSeqKey[fastaSeqKeyShort.index(seqID)]])
                    pbar.update()
                if len(localFailedDownload):
                    logger.warning(f"Unable to download these sequences: {', '.join(localFailedDownload)}")
    
    if len(failedDownload):
        logger.info("Summary of undownloaded IDs")
        displaySeqIds(failedDownload)
    else:
        logger.info("The download went well")

