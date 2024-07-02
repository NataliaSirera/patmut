import tempfile
from collections import OrderedDict
from Bio.SeqIO import parse, write
from io import StringIO


from config import shell


def searchHomologous(blastAlgPath, blastProtDBPath, numIter, evalue, protSeq):
    numCores, error = shell("nproc", errDiscard=False)
    if error:
        numCores, error = shell("sysctl -n hw.ncpu", errDiscard=False)
        if error:
            numCores = 4
    tmpFile = tempfile.NamedTemporaryFile()
    with open(tmpFile.name, "w") as fw:
        fw.write(">id\n" + protSeq)
    blastSearch = shell(
        blastAlgPath
        + "psiblast -db "
        + blastProtDBPath
        + " -query "
        + tmpFile.name
        + " -num_iterations "
        + str(numIter)
        + " -evalue "
        + str(evalue)
        + ' -outfmt "7 sacc slen nident" -num_threads '
        + numCores
    )

    flag = False
    lastRoundBlastSearch = OrderedDict()
    for line in blastSearch.split("\n"):
        if "Iteration: " + str(numIter) in line:
            flag = True
        elif flag:
            if line.startswith("#") or line == "" or line.startswith("Search has"):
                continue
            line = line.split("\t")
            seqID = line[0]
            try:
                identities = [lastRoundBlastSearch[seqID][1], float(line[2])]
                iden = max(identities)
                index = identities.index(iden)
                length = [lastRoundBlastSearch[seqID][0], int(line[1])][index]
                lastRoundBlastSearch[seqID] = [length, iden]
            except KeyError:
                lastRoundBlastSearch[seqID] = [int(line[1]), float(line[2])]
    return lastRoundBlastSearch


def filterHomologous(foundHomologous, qLen, minSeqIden):
    filteredHomologous = [
        sID
        for sID, [sLen, iden] in foundHomologous.items()
        if iden / ((qLen + sLen) / 2.0) * 100 >= float(minSeqIden)
    ]
    return filteredHomologous


def retrieveSeqFromID(protDBPath, idList, uniprotID, blastAlgPath, protSeq):
    tmpFile = tempfile.NamedTemporaryFile()
    with open(tmpFile.name, "w") as fw:
        fw.write("\n".join(idList))
    seqDBString = shell(
        blastAlgPath + "blastdbcmd -db " + protDBPath + " -entry_batch " + tmpFile.name
    )
    tmpFile.close()

    # Check MSA has the query protein and it has not been lost during MSA building (derive)
    if "UniRef100_" + uniprotID not in idList:
        seqDBString += "\n" + ">lcl|UniRef100_" + uniprotID + "\n" + protSeq

    seqDB = list(parse(StringIO(seqDBString), "fasta"))

    return seqDB


def buildMsaMuscle(
    filteredHomoProtID, muscleAlgPath, blastProtDBPath, uniprotID, blastAlgPath, protSeq
):
    homologousDB = retrieveSeqFromID(
        blastProtDBPath, filteredHomoProtID, uniprotID, blastAlgPath, protSeq
    )
    tmpFile = tempfile.NamedTemporaryFile()
    with open(tmpFile.name, "w") as fw:
        write(homologousDB, fw, "fasta")

    homologousAlignment = shell(muscleAlgPath + "muscle -in " + tmpFile.name)
    msaDB = list(parse(StringIO(homologousAlignment), "fasta"))
    tmpFile.close()

    return msaDB


def buildMSA(
    protSeq, blastAlgPath, blastProtDBPath, numIter, evalue, minSeqIden, muscleAlgPath, uniprotID
):
    foundHomologous = searchHomologous(blastAlgPath, blastProtDBPath, numIter, evalue, protSeq)
    filteredHomologous = filterHomologous(foundHomologous, len(protSeq), minSeqIden)
    msaDB = buildMsaMuscle(
        filteredHomologous, muscleAlgPath, blastProtDBPath, uniprotID, blastAlgPath, protSeq
    )

    return msaDB
