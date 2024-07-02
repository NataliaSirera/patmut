from math import log
from Bio.Data.IUPACData import extended_protein_letters


def henikoffWeights(msaDB, msaDBT):
    width, length = len(msaDBT), len(msaDBT[0])
    colAA = [
        dict((aa, col.count(aa)) for aa in list(extended_protein_letters) + ["-"])
        for pos, col in enumerate(msaDBT)
    ]
    nonZeros = [len([aa for aa, count in col.items() if count > 0]) for col in colAA]
    colWeights = [
        dict(
            (aa, 1.0 / (colAA[pos][aa] * nonZeros[pos])) if colAA[pos][aa] > 0 else (aa, 0.0)
            for aa in list(extended_protein_letters) + ["-"]
        )
        for pos in range(len(colAA))
    ]
    recWeights = [
        (sum([colWeights[pos][aa] for pos, aa in enumerate(str(record.seq))])) / float(width)
        for record in msaDB
    ]
    return recWeights


def ShannonEntropy(msaDB, msaDBT, alignedQuery, weights, variants):
    colEntropy = [
        dict(
            (aa, sum([weights[seq] for seq, msaAA in enumerate(col) if msaAA == aa]))
            for aa in list(extended_protein_letters) + ["-"]
        )
        for pos, col in enumerate(msaDBT)
        if alignedQuery[pos] != "-"
    ]
    recEntropy = [
        -sum([aaE * log(aaE, 2) for aaE in col.values() if aaE > 0]) for col in colEntropy
    ]
    entropyARFF = dict((variant, recEntropy[int(variant[1:-1]) - 1]) for variant in variants)
    return entropyARFF, recEntropy


def entropy(msaDB, msaDBT, alignedQuery, variants):
    weights = henikoffWeights(msaDB, msaDBT)
    entropyARFF, recEntropy = ShannonEntropy(msaDB, msaDBT, alignedQuery, weights, variants)
    return entropyARFF, recEntropy


def entropyAvgWindow(recEntropy, variants):
    # entropy average window
    newEntropy = []
    left_window = []
    right_window = []

    for index in range(len(recEntropy)):
        if index < 3:  # for cases where less than 3 residue on the left
            left_window = recEntropy[:index]
            right_window = recEntropy[index + 1 : index + 4]
            newEntropy.append(
                (sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window))
            )
        elif len(recEntropy) - index < 3:  # for cases where less than 3 residue on the right
            left_window = recEntropy[index - 3 : index]
            right_window = recEntropy[index + 1 : len(recEntropy)]
            newEntropy.append(
                (sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window))
            )
        else:
            left_window = recEntropy[index - 3 : index]
            right_window = recEntropy[index + 1 : index + 4]
            newEntropy.append(
                (sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window))
            )

    entropyAvgWindow = dict((variant, newEntropy[int(variant[1:-1]) - 1]) for variant in variants)

    return entropyAvgWindow, newEntropy
