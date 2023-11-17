import pandas, os
from alive_progress import alive_bar
from Bio import Entrez

def searchIndex(idList, identifier, keyLenght, indexStart, IndexEnd):
    startString = identifier[:keyLenght]
    endString = identifier[-keyLenght:]
    for key in indexStart[startString]:
        if key in IndexEnd[endString] and identifier==idList[key]:
            return key
    return None

def ParseSequences(stream):
    print("Parsing data")
    sequences = {}
    id = None
    buffer = []
    for line in stream:
        if ">" in line:
            if id:
                sequences[id] = buffer
                buffer = []
            id = line[1:].split(" ")[0]
            buffer.append(line)
        else:
            buffer.append(line)
    sequences[id] = buffer
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

def searchEntrez(query):
    print("Searching Entrez IDs")
    subQueryString = [i[1] for i in query]
    queryString = ",".join(subQueryString)
    # print(queryString)
    # handle = Entrez.esearch(db="nuccore", rettype="Fasta", retmode="text", term=queryString)
    handle = Entrez.esearch(db="nuccore", retmode="text", term=queryString)
    records = Entrez.read(handle)
    handle.close()
    return records["IdList"]

def fetchEntrez(entrezIDs):
    print("Fetching sequences")
    handle = Entrez.efetch(db="nuccore", id=entrezIDs, rettype="fasta", retmode="text")
    data = handle.readlines()
    handle.close()
    return data



indexFile = "/mnt/d/bioDB/ProMGE/mges.txt"
dbLocation = "/mnt/d/bioDB/ProMGE/sequencesData"
Entrez.email="antoine.druart@me.com"
# indexSequenceFile = "/mnt/h/bioDB/ProMGE/indexSequences.txt"

print(f"Reading {indexFile}")
data = pandas.read_table(indexFile, sep="\t")

if not os.path.exists(dbLocation):
    os.mkdir(dbLocation)

seqIDs = []
count = 0
indexSize = 3
startIndex = {}
endIndex = {}
query = []
failedDownload = []
numberOfSeqPerDl = 21
nbseq = len(data["mge_genome_position"])
with alive_bar(nbseq) as bar: 
    for i in range(nbseq):
        seqID = data["mge_genome_position"][i].split(":")[0].split(".")[2]
        projectID = data["mge_genome_position"][i].split(":")[0].split(".")[1]
        startString = seqID[:indexSize]
        endString = seqID[-indexSize:]
        if startString not in startIndex.keys():
            startIndex[startString] = []
        if endString not in endIndex.keys():
            endIndex[endString] = []
        # print(seqIDs, seqID, indexSize, startIndex, endIndex)
        itemIndex = searchIndex(seqIDs, seqID, indexSize, startIndex, endIndex)
        if itemIndex is None and seqID not in seqIDs:
            seqIDs.append(seqID)
            startIndex[startString].append(len(seqIDs)-1)
            endIndex[endString].append(len(seqIDs)-1)
            query.append(seqID)
        bar()
        # print(len(seqIDs))
        # if seqID not in query:
        if len(query)>=numberOfSeqPerDl or i==nbseq-1 and len(query)>0:
            entrezIDs = searchEntrez(query)
            print(f"EntrezIDs : {len(entrezIDs)}")
            rawSequences = fetchEntrez(entrezIDs)
            print(f"rawSequencesIDs : {rawSequences}")
            FastaSequences = ParseSequences(rawSequences)
            print(f"EntrezIDs : {len(FastaSequences.keys())}")
            FastaSequencesKeysShort = []
            FastaSequencesKeysLong = [seqName for seqName in FastaSequences.keys()]
            for sID in query:
                # print(sID)
                # print(FastaSequences.keys())
                if sID[1] in FastaSequencesKeysShort:
                    print(f"Saving {sID[1]}")
                    seqPath = getDestinationFolder(sID[1], dbLocation)
                    seqIndex = FastaSequencesKeysLong[FastaSequencesKeysShort.index(sID[1])]
                    writeSequence(seqPath, sID, FastaSequences[seqIndex])
                    print(f"{sID[1]} saved")
                else:
                    failedDownload.append(sID)
                query = []
            if len(failedDownload):
                print(f"\033[31mFailed to download ({len(failedDownload)}): {failedDownload}\033[0m")
            failedDownload = []
print("Fin du téléchargement")