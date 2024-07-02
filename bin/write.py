import os
import Bio.SeqIO

from config import version


def roundAttribute(value):
    try:
        return "%.3f" % round(float(value), 3)
    except ValueError:
        return value


def writeMsa(msaDB, outPath):
    with open(outPath + ".msa", "w") as fw:
        Bio.SeqIO.write(msaDB, fw, "fasta")


def writeNeutresMsaID(neutralMSAID, outPath):
    with open(outPath + ".neutresMSA", "w") as fw:
        fw.write(
            "\n".join(
                [
                    "\t".join([variant, seqID])
                    for variant, seqID in sorted(
                        neutralMSAID, key=lambda x: (int(x[0][1:-1]), x[0][-1])
                    )
                ]
            )
        )


def writeArff(configOptions, outPath):
    features = [
        "vdwVolume",
        "hydrophobicity",
        "substitutionMatrix",
        "pssm-native",
        "pssm-mutated",
        "entropy",
        "freqGaps",
        "entropyAvgWindow",
        "relEntropy",
        "fractionSimilarity",
        "impRes",
    ]
    fTypes = [
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "NUMERIC",
        "{True,False}",
    ]
    userFeatures = [f for f in features if configOptions[f].value]

    arff = "@RELATION " + os.path.abspath(configOptions["config"].value) + "_v" + version() + "\n\n"
    arff += "@ATTRIBUTE variant STRING\n@ATTRIBUTE uniprotID STRING\n"
    arff += "".join(
        [
            "@ATTRIBUTE " + feature + " " + fType + "\n"
            for feature, fType in zip(features, fTypes)
            if configOptions[feature].value
        ]
    )
    arff += "@ATTRIBUTE tag {1,0}\n"
    arff += (
        "\n@DATA\n"
        + "\n".join(
            ",".join(
                [variant, configOptions["uniprotID"].value]
                + [
                    roundAttribute(configOptions[f].value[variant])
                    for f in userFeatures
                    if configOptions[f].value
                ]
                + [configOptions["tag"].value[variant]]
            )
            for variant in configOptions["variants"].value
        )
        + "\n"
    )
    with open(outPath + ".arff", "w") as fw:
        fw.write(arff)


def write(configOptions):
    outPath = configOptions["output"].value
    if "msa" in configOptions["write"].value:
        writeMsa(configOptions["msaDB"].value, outPath)
    if "neutresMsaID" in configOptions["write"].value:
        writeNeutresMsaID(configOptions["neutralMSAID"].value, outPath)
    if "arff" in configOptions["write"].value:
        writeArff(configOptions, outPath)
