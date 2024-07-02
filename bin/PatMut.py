import config
import impRes
import vdwVolume
import hydrophobicity
import substitutionMatrix
import msa
import pssm
import entropy
import freqGaps
import write


if __name__ == "__main__":

    # Set Up
    configOptions = config.config()

    # Build MSA
    if configOptions["msaBuild"].value:
        configOptions["msaDB"].value = msa.buildMSA(
            configOptions["protSeq"].value,
            configOptions["blastPath"].value,
            configOptions["blastProtDB"].value,
            configOptions["numIter"].value,
            configOptions["evalue"].value,
            configOptions["minSeqIden"].value,
            configOptions["muscleAlg"].value,
            configOptions["uniprotID"].value,
        )
        if "arff" in configOptions["write"].user or "neutresMsaID" in configOptions["write"].user:
            # print(configOptions["msaDB"].value)
            # print(configOptions["uniprotID"].value)
            configOptions["alignedQuery"].value = config.retrieveAlignedQuery(
                configOptions["msaDB"].value,
                configOptions["uniprotID"].value,
                configOptions["protSeq"].value,
            )

    # Parse MSA variants
    if "neutresMsaID" in configOptions["write"].user or (
        configOptions["neutral"].user
        and (
            "msaDB" in configOptions["neutral"].user or "msaBuild" in configOptions["neutral"].user
        )
    ):
        configOptions = config.retrieveMSANeutralVariants(configOptions)

    # Retrieve all variants
    variants = sorted(
        [
            (variant, tag)
            for varPop, tag in zip(["pathological", "neutral", "predicting"], ["1", "0", "?"])
            if configOptions[varPop].user
            for variant in configOptions[varPop].value
        ],
        key=lambda x: (int(x[0][1:-1]), x[0][-1]),
    )
    configOptions["variants"].value = [i[0] for i in variants]
    configOptions["tag"].value = dict((i[0], i[1]) for i in variants)

    # Parse Sequence Features
    if configOptions["impRes"].value:
        configOptions["impRes"].value = impRes.impRes(
            configOptions["impResDB"].value, configOptions["variants"].value
        )
    if configOptions["vdwVolume"].value:
        configOptions["vdwVolume"].value = vdwVolume.vdwVolume(
            configOptions["vdwVolumeDB"].value,
            configOptions["variants"].value,
            configOptions["absoluteFlag"].value,
        )
    if configOptions["hydrophobicity"].value:
        configOptions["hydrophobicity"].value = hydrophobicity.hydrophobicity(
            configOptions["hydrophobicityDB"].value,
            configOptions["variants"].value,
            configOptions["absoluteFlag"].value,
        )
    if configOptions["substitutionMatrix"].value:
        configOptions["substitutionMatrix"].value = substitutionMatrix.substitutionMatrix(
            configOptions["substitutionMatrixDB"].value, configOptions["variants"].value
        )

    # Parse MSA Features
    if any(
        [
            configOptions["pssm-native"].value,
            configOptions["pssm-mutated"].value,
            configOptions["entropy"].value,
            configOptions["freqGaps"].value,
        ]
    ):
        configOptions["msaDBT"].value = [
            "".join(col)
            for col in zip(*[str(record.seq) for record in configOptions["msaDB"].value])
        ]
    if any([configOptions["pssm-native"].value, configOptions["pssm-mutated"].value]):
        pssmValue, RE_all = pssm.pssm(
            configOptions["msaDBT"].value,
            configOptions["alignedQuery"].value,
            configOptions["protSeq"].value,
        )
    if configOptions["pssm-native"].value:
        configOptions["pssm-native"].value = pssm.pssmn(pssmValue, configOptions["variants"].value)
    if configOptions["pssm-mutated"].value:
        configOptions["pssm-mutated"].value = pssm.pssmm(pssmValue, configOptions["variants"].value)
    if configOptions["entropy"].value:
        configOptions["entropy"].value, recEntropy = entropy.entropy(
            configOptions["msaDB"].value,
            configOptions["msaDBT"].value,
            configOptions["alignedQuery"].value,
            configOptions["variants"].value,
        )
    if configOptions["freqGaps"].value:
        configOptions["freqGaps"].value = freqGaps.freqGaps(
            configOptions["msaDBT"].value,
            configOptions["alignedQuery"].value,
            configOptions["variants"].value,
        )

    if configOptions["entropyAvgWindow"].value:
        configOptions["entropyAvgWindow"].value, newEntropy = entropy.entropyAvgWindow(
            recEntropy, configOptions["variants"].value
        )
    if configOptions["relEntropy"].value:
        configOptions["relEntropy"].value = pssm.relativeEntropy(
            RE_all, configOptions["variants"].value
        )
    if configOptions["fractionSimilarity"].value:
        configOptions["fractionSimilarity"].value = pssm.fractionSimilarity(
            configOptions["msaDB"].value,
            configOptions["alignedQuery"].value,
            configOptions["variants"].value,
        )

    # Write Output
    write.write(configOptions)
