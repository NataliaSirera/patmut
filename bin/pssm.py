from math import log
from Bio.Data.IUPACData import protein_letters
from collections import OrderedDict
import numpy as np


def pssm(msaDBT, alignedQuery, protSeq):
    width, length, aaDB = len(msaDBT), len(msaDBT[0]), "".join([col for col in msaDBT])
    colAA = [
        dict((aa, col.count(aa) / float(length)) for aa in list(protein_letters))
        for pos, col in enumerate(msaDBT)
    ]
    totalAA = dict((aa, aaDB.count(aa) / float(length * width)) for aa in list(protein_letters))
    pssm = [
        dict(
            (
                aa,
                log((colAA[pos][aa]) / (totalAA[aa]), 2)
                if count != 0
                else log(width / (length * width - (length - 1.0)), 2),
            )
            for aa, count in colAA[pos].items()
        )
        for pos, col in enumerate(msaDBT)
        if alignedQuery.seq[pos] != "-"
    ]
    rpssm = [
        dict(
            (
                aa,
                (colAA[pos][aa]) * log((colAA[pos][aa]) / (totalAA[aa]), 2)
                if count != 0
                else 0 * log(width / (length * width - (length - 1.0)), 2),
            )
            for aa, count in colAA[pos].items()
        )
        for pos, col in enumerate(msaDBT)
        if alignedQuery.seq[pos] != "-"
    ]

    RE_all = []
    left_window = []
    right_window = []

    for index in range(len(rpssm)):
        aa = protSeq[index]
        if index < 3:
            left_window = [pos[aa] for pos in rpssm[:index]]
            right_window = [pos[aa] for pos in rpssm[index + 1 : index + 4]]
            RE_all.append(
                ((sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window)))
                + rpssm[index][aa]
            )
        elif len(rpssm) - index < 3:
            left_window = [pos[aa] for pos in rpssm[index - 3 : index]]
            right_window = [pos[aa] for pos in rpssm[index + 1 : len(rpssm)]]
            RE_all.append(
                ((sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window)))
                + rpssm[index][aa]
            )
        else:
            left_window = [pos[aa] for pos in rpssm[index - 3 : index]]
            right_window = [pos[aa] for pos in rpssm[index + 1 : index + 4]]
            RE_all.append(
                ((sum(left_window) + sum(right_window)) / (len(left_window) + len(right_window)))
                + rpssm[index][aa]
            )

    return pssm, RE_all


def pssmn(pssm, variants):
    pssmnARFF = dict((variant, pssm[int(variant[1:-1]) - 1][variant[0]]) for variant in variants)
    return pssmnARFF


def pssmm(pssm, variants):
    pssmmARFF = dict((variant, pssm[int(variant[1:-1]) - 1][variant[-1]]) for variant in variants)
    return pssmmARFF


def relativeEntropy(RE_all, variants):
    relEntropy = dict((variant, RE_all[int(variant[1:-1]) - 1]) for variant in variants)
    return relEntropy


def fractionSimilarity(msaDB, alignedQuery, variants):
    fracSim = []

    fractionDict = dict()
    for seq in msaDB:
        any_pair = 0
        exact_pair = 0
        for pos, human_res in enumerate(alignedQuery):
            if human_res != "-":
                other_res = seq[pos]
                if other_res != "-":
                    any_pair += 1
                if human_res == other_res:
                    exact_pair += 1
        fractionDict[str(seq.seq)] = exact_pair / any_pair
    fractionDict_sorted = OrderedDict(
        (sorted(fractionDict.items(), key=lambda item: item[1], reverse=True))
    )

    median_all = np.median([fractionDict_sorted[k] for k in list(fractionDict_sorted)])

    for pos, human_res in enumerate(alignedQuery):
        first_n = []
        exact_n = []
        n = 0
        if human_res != "-":
            for seq, frac in fractionDict_sorted.items():
                other_res = seq[pos]
                if human_res == other_res:
                    exact_n.append(frac)
                    n += 1
            first_n = [fractionDict_sorted[k] for k in list(fractionDict_sorted)[:n]]

            if n == 0:
                fracSim.append(0)
            else:
                fracSim.append(abs(np.median(first_n) - np.median(exact_n)) / median_all)
    fractionSimilarities = dict((variant, fracSim[int(variant[1:-1]) - 1]) for variant in variants)
    return fractionSimilarities
