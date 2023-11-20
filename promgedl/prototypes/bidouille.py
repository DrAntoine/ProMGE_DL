import os, re
# from alive_progress import alive_bar
from tqdm import tqdm
from Bio import Entrez

class noeud:

    def __init__(self, IdPere, value, idFilsG=None, idFrereD=None) -> None:
        self.pere = IdPere
        self.content=value
        self.filsG = idFilsG
        self.frereD = idFrereD

def makeFilteredTree(sequences):
    tree = [noeud(IdPere=None, value="START")]
    nbseq = len(sequences)
    # with alive_bar(nbseq) as bar:
    with tqdm(total=nbseq, desc="Sorting accession numbers (1/2)") as pbar:
        for i in range(nbseq):
            # seqID = data["mge_genome_position"][i].split(":")[0].split(".")[2]
            seqID = sequences[i]
            # seqID+="!"
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
                    while childNode.content != char :# and childNode.frereD != None :
                        if childNode.frereD is None:
                            tree.append(noeud(IdPere=nodeIndex, value=char))
                            childNode.frereD = len(tree)-1
                        else:
                            childNodeindex = childNode.frereD
                        childNode = tree[childNodeindex]
                    nodeIndex=childNodeindex
            # bar()
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
                # if node.content != "!":
                seq+=node.content
                node = tree[node.pere]
            seq = seq[::-1]
            # print(seq)
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

def ParseSequences(stream):
    print("Parsing data")
    sequences = {}
    seqid = None
    buffer = []
    print(stream)
    input("wait")
    for line in stream:
        if ">" in line:
            if id:
                sequences[id] = buffer
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

def searchEntrez(query, entrezDB = "nuccore"):
    print(f"Searching Entrez IDs db:{entrezDB}")
    subQueryString = list(query)
    queryString = ", ".join(subQueryString)
    # print(queryString)
    # handle = Entrez.esearch(db="nuccore", rettype="Fasta", retmode="text", term=queryString)
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

def __count_generator__(reader):
    b = reader(1024*1024)
    while b:
        yield b
        b = reader(1024*1024)

def readFile(filePath, filters):
    print(f"Reading {filePath}")
    sequencesID = []
    dbindex =[]
    lineCount = 0
    print("Estimating fileLenght")
    with open(filePath, "rb") as dbfile:
        count_generator = __count_generator__(dbfile.raw.read)
        lineCount = sum(buffer.count(b'\n') for buffer in count_generator)
    # print(lineCount)
    currentLine = 0
    # with alive_bar(lineCount) as bar:
    with tqdm(total=lineCount, desc="Extracting accession numbers") as pbar:
        with open(filePath, "r") as dbfile:
            for line in dbfile:
                lineSplited = line.split("\t")
                if currentLine == 0:
                    dbindex = lineSplited
                else:
                    acceptLine = True
                    for filter in filters.keys():
                        if filters[filter] != ["*"] and lineSplited[dbindex.index(filter)] not in filters[filter]:
                            acceptLine = False
                            break
                    if acceptLine:
                        sequencesID.append(lineSplited[0].split(":")[0].split(".")[2])
                # bar()
                pbar.update()
                currentLine+=1
    return sequencesID

def downloadSeq(seqIDs, dbLocation, entrezDB = "nucleotide"):
    fullListOfFailed = []
    entrezDB = "nuccore"
    # with alive_bar(len(seqIDs)) as bar:
    with tqdm(total=len(seqIDs), desc="Downloading sequences") as pbar:
        query = []
        numberOfSeqPerDl = 20
        failedDownload = []

        for i in range(len(seqIDs)): #seqID in seqIDs:
            query.append(seqIDs[i])
        # print(len(seqIDs))
            # if seqID not in query:
            if len(query)>=numberOfSeqPerDl or i==len(seqIDs)-1 and len(query)>0:
                succes = False
                FastaSequences = None

                for trydl in range(3):
                    try:
                        entrezIDs = searchEntrez(query, entrezDB)
                        if entrezIDs == []: raise ValueError("Idlist problem")
                        rawSequences = fetchEntrez(entrezIDs, entrezDB)
                        FastaSequences = ParseSequences(rawSequences)
                    except RuntimeError as error:
                        print("Un problème est survenu !")
                        print(error)
                    except:
                        print("Un problème est survenu !")
                        print("Erreur non gérée")
                    if FastaSequences is not None:
                        succes = True
                        break
                    else:
                        print(f'Try {trydl+1} failed')
                
                print(FastaSequences)
                if not succes:
                    for sID in query:
                        failedDownload.append(sID)
                        # bar()
                        pbar.update()
                else:
                    for seqName in FastaSequences.keys():
                        print(seqName)
                    FastaSequencesKeysShort = [seqName.split(".")[0] for seqName in FastaSequences.keys() if seqName is not None ]
                    FastaSequencesKeys = list(FastaSequences.keys())
                # print(FastaSequencesKeys)
                    for sID in query:
                        # print(sID)
                        # print(FastaSequences.keys())
                        if sID in FastaSequencesKeysShort:
                            # print(f"Saving {sID}")
                            seqPath = getDestinationFolder(sID, dbLocation)
                            # seqIndex = FastaSequencesKeysLong[FastaSequencesKeysShort.index(sID[1])]
                            writeSequence(seqPath, sID, FastaSequences[FastaSequencesKeys[FastaSequencesKeysShort.index(sID)]])
                            # print(f"{sID} saved")
                        else:
                            failedDownload.append(sID)
                        # bar()
                        pbar.update()
                string = "ID :"
                for q in query:
                    if q in failedDownload:
                        string += f" \033[31m{q}\033[0m,"
                    else:
                        string += f" \033[92m{q}\033[0m,"
                print(string[:-1])
                if len(failedDownload):
                    fullListOfFailed+=failedDownload
                failedDownload = []
                query = []
    return fullListOfFailed


indexFile = "/mnt/d/bioDB/ProMGE/mges.txt"
# indexFile = "H:/bioDB/ProMGE/mges.txt"
dbLocation = "/mnt/d/bioDB/ProMGE/sequencesData"
# dbLocation = "H:/bioDB/ProMGE/sequencesData"

# indexSequenceFile = "/mnt/h/bioDB/ProMGE/indexSequences.txt"

# data = pandas.read_table(indexFile, sep="\t")
filters = {
    "mge_genome_position":["*"],
    "genus": ["Salmonella"]
}
sequences = readFile(indexFile, filters)
if len(sequences)== 0:
    print("The filters are to strict. No sequence availble")
    exit()

# if not os.path.exists(dbLocation):
    # os.mkdir(dbLocation)

sequencesID = extractSequenceIDFromTree(makeFilteredTree(sequences))
sequencesID.sort()

# print(sequencesID)
# print(len(sequencesID))
# step = 10
# for i in range((len(sequencesID)//step)+1):
#     print(sequencesID[i*step:(i+1)*step])
GBNucleotids1 = re.compile("^[A-Z]{1}[0-9]{5}")
GBNucleotids2 = re.compile("^[A-Z]{2}[0-9]{6}")
GBNucleotids3 = re.compile("^[A-Z]{2}[0-9]{8}")
GBProteins1 = re.compile("^[A-Z]{3}[0-9]{5}")
GBProteins2 = re.compile("^[A-Z]{3}[0-9]{7}")
GBwgs1 = re.compile("^[A-Z]{4}[0-9]{2}[0-9]{6}[0-9]*")
GBwgs2 = re.compile("^[A-Z]{6}[0-9]{2}[0-9]{7}[0-9]*")
GBMGA = re.compile("^[A-Z]{5}[0-9]{7}")

# refseq1 = re.compile("^AC_.+")

GB_nuc = []
GB_prot = []
GB_MGA = []
GB_WGS = []
other = []

for seq in sequencesID:
    if GBNucleotids1.match(seq) or GBNucleotids2.match(seq) or GBNucleotids3.match(seq):
        GB_nuc.append(seq)
    elif GBProteins1.match(seq) or GBProteins2.match(seq):
        GB_prot.append(seq)
    elif GBMGA.match(seq):
        GB_MGA.append(seq)
    elif GBwgs1.match(seq) or GBwgs2.match(seq):
        GB_WGS.append(seq)
    else:
        other.append(seq)

print(f"GB_nuc : {len(GB_nuc)}")
print(f"GB_prot : {len(GB_prot)}")
print(f"GB_MGA : {len(GB_MGA)}")
print(f"GB_WGS : {len(GB_WGS)}")
print(f"OTHER : {len(other)}")

Failed_GB_nuc = downloadSeq(GB_nuc, dbLocation)
Failed_GB_prot = downloadSeq(GB_prot, dbLocation)
Failed_GB_MGA = downloadSeq(GB_MGA, dbLocation)
Failed_GB_WGS = downloadSeq(GB_WGS, dbLocation)
Failed_other = downloadSeq(other, dbLocation)

print(f"Failed_GB_nuc : {len(Failed_GB_nuc)}/{len(GB_nuc)}")
print(f"Failed_GB_prot : {len(Failed_GB_prot)}/{len(GB_prot)}")
print(f"Failed_GB_MGA : {len(Failed_GB_MGA)}/{len(GB_MGA)}")
print(f"Failed_GB_WGS : {len(Failed_GB_WGS)}/{len(GB_WGS)}")
print(f"Failed_OTHER : {len(Failed_other)}/{len(other)}")


print("Fin du téléchargement")