import re
import os
import sys
import subprocess
from Bio.Data.IUPACData import (
    protein_letters,
    protein_letters_3to1,
    extended_protein_letters,
)
from Bio.Align import substitution_matrices
import Bio.SeqIO
from urllib.request import urlopen
from urllib.error import HTTPError


class option:
    def __init__(self, dataType=None, options=None, dependencies=None, default=None, user=None):

        """
        Class: option
        Attributes:
            - Type: can be boolean, regex, path, float or list
            - Options: values that user can introduce for this config option
            - Dependencies: options needed to compute this option
            - User: user value of the option
            - Value: final value of the option
        """

        self.type = dataType
        self.options = options
        self.dependencies = dependencies
        self.user = user
        self.value = default


def shell(cmd, errDiscard=True):
    """
    Input: command line to execute
    Defaults Param: errDiscard=True
    Function: command line is executed in shell, stdout and stderr are redirected
    Output:
            - Default mode: stdout only
            - errDiscard=False: stdout,stderr.rstrip('\n')
    """

    stdout, stderr = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    stdout, stderr = stdout.decode("ascii").rstrip("\n"), stderr.decode("ascii").rstrip("\n")
    if errDiscard:
        return stdout
    else:
        return stdout, stderr


def file2dict(filename):
    """
    Input: file path
    Format: key \t value \n
    Function: retrieve data from file to dictionary
    Output: dict
    """

    data = dict()
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip("\n").split("\t")
            data[line[0]] = line[1]
    return data


def checkRes(residue, protLength):
    """
    Input: residue and protein length
    Format: 1 based coordinates
    Function: check if the residue is between 1 and the protein length
    Output: True/False
    """
    try:
        residue = int(residue)
    except ValueError:
        return False
    if 1 <= residue <= protLength:
        return True
    else:
        return False


def checkNaturalAa(seq, extended=False):
    """
    Input: string of amino acid(s)
    Format: one letter based
    Function: check if the amino acid is natural or not
    Output: True/False
    """

    naturalAa = "^(" + "|".join(list(protein_letters)) + ")+$"
    if re.search(naturalAa, seq):
        return True
    else:
        return False


def checkProtAa(variant, protSeq):
    """
    Input: variant to check and protein sequence
    Format: one letter based both variant and protein, variant coordinates are 1-based
    Function: check if both variant's original and mutated amino acids are natural amino acids,
              if the variant's residue is within the protein and if the variant's original amino acid
              corresponds to the protein's amino acid at the variant's residue position
    Output: True/False
    """

    oriAa = variant[0]
    residue = variant[1:-1]
    mutAa = variant[-1]
    if checkNaturalAa(oriAa) and checkNaturalAa(mutAa) and checkRes(residue, len(protSeq)):
        if protSeq[int(residue) - 1] == oriAa:
            return True
    return False


def cleanVariant(variant, protSeq):
    """
    Input: variant
    Format: residue: 1 based coordinates, amino acids: 1 or 3 letters
    Function: convert amino acid from 3 to 1 and check if amino acids of variant are natural amino acids
              and corresponds to protSeq and residue is within protein length
    Output: clean variant / None
    """

    variantFormat = re.search("([a-zA-Z]+)(\d+)([a-zA-Z]+)", variant)
    if not variantFormat or len(variantFormat.group(1)) not in [1, 3]:
        return None
    if len(variantFormat.group(1)) == 3:
        try:
            oriAa = protein_letters_3to1[variant[:3]]
        except KeyError:
            return None
        try:
            mutAa = protein_letters_3to1[variant[-3:]]
        except KeyError:
            return None
        variant = oriAa + variant[3:-3] + mutAa
    if checkProtAa(variant, protSeq):
        return variant
    else:
        return None


def retrieveProtSeq(uniprotID, blastProtDB, blastPath):
    """
    Input: uniprot accession number, blast database of uniprot proteins, path to blast software
    Function: blastdbcmd to find the protein sequence of the uniprotID in the protein blast database
    Output: protSeq of uniprotID / sys.exit() if stderr in blastdbcmd
    """

    stdout, stderr = shell(
        blastPath + "blastdbcmd -db " + blastProtDB + " -entry UniRef100_" + uniprotID,
        errDiscard=False,
    )
    if stderr:
        if "not found" in stderr:
            try:
                html = (
                    urlopen("http://www.uniprot.org/uniprot/" + uniprotID + ".fasta")
                    .read()
                    .decode("utf-8")
                )
            except HTTPError:
                sys.exit(
                    "ERROR: "
                    + uniprotID
                    + " not found neither blastDB nor http://www.uniprot.org/uniprot/"
                    + uniprotID
                )
            else:
                if html == "":
                    sys.exit(
                        "ERROR: seems to be a problem in uniprot, please retrieve the sequence yourself."
                        + uniprotID.split("-")[0]
                    )
                else:
                    return "".join(html.split("\n")[1:])
        else:
            sys.exit("ERROR: blastdbcmd " + stderr)
    else:
        return "".join(stdout.split("\n")[1:])


def retrieveImpRes(filename, uniprotID, protLength):
    """
    Input: file path, uniprotID/'', protein Length
    Function: retrieve important residues from humsavar (uniprotID given, column=1) or
              userImpResDB (uniprotID='',column=0). Check if residue is between 1 and protein length
    Output: list of important residues / Discarded if check not passed
    """

    impResVars = []
    if uniprotID:
        with open(filename) as fh:
            for line in fh:
                if line.startswith(uniprotID):
                    residues = line.rstrip("\n").split("\t")[1].split(",")
                    for residue in residues:
                        if checkRes(residue, protLength):
                            impResVars.append(residue)
                        else:
                            print(
                                "WARNING: important residue "
                                + str(residue)
                                + " discarded because it is not between 1 and "
                                + str(protLength)
                            )
                    break
    else:
        with open(filename) as fh:
            for residue in fh.read().rstrip("\n").split("\n"):
                if checkRes(residue, protLength):
                    impResVars.append(residue)
                else:
                    print(
                        "WARNING: important residue "
                        + str(residue)
                        + " discarded because it is not between 1 and "
                        + str(protLength)
                    )
    return impResVars


def retrieveUserVariants(filename, protSeq):
    """
    Input: file path with a variant per line
    Format: 3 or 1 letter code is possible for amino acids, residue are 1 based
    Function: retrieve variants from a file, one per line, and convert it to 1 letter code if need it.
              Check natural amino acid, original amino acid protein and residue protein length
    Output: list of variants / Discarded if check not passed
    """

    userVars = []
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                variant = cleanVariant(line, protSeq)
                if variant:
                    userVars.append(variant)
                else:
                    print(
                        "WARNING: variant "
                        + line
                        + " discarded because it does not match protein sequence"
                    )
    return userVars


def retrieveHumsavarVariants(filename, uniprotID, protSeq):
    """
    Input: file path, uniprotID, protein sequence
    Function: retrieve disease and polymorphism variants as pathological and neutral variants from humsavar database.
              Check if natural amino acid, original amion acid protein and residue protein length.
    Output: pathological variants, neutral variants from humsavar
    """

    pathoVariants = []
    neuVariants = []
    with open(filename) as fh:
        for line in fh:
            if uniprotID in line:
                line = line.split()
                variant = cleanVariant(line[3].lstrip("p."), protSeq)
                if variant:
                    if line[4] == "Disease":
                        pathoVariants.append(variant)
                    elif line[4] == "Polymorphism":
                        neuVariants.append(variant)
                else:
                    pass
    return pathoVariants, neuVariants


def retrieveSysMutVariants(protSeq):
    """
    Input: protein sequence
    Function: calculate all possible non synonymous variants for each residue of the protein sequence
    Output: list of all non synonymous variants for a protein
    """

    sysMuts = []
    for ind, oriAA in enumerate(list(protSeq)):
        for mutAA in protein_letters:
            if oriAA != mutAA:
                sysMuts.append(oriAA + str(ind + 1) + mutAA)
    return sysMuts


def pathoNeuPredPriority(configOptions):
    """
    Input: config options
    Function: remove variants that are in another set, with priority Pathological > Neutral > Predicting
    Output: config options
    """

    for variantPopulation in ["pathological", "neutral", "predicting"]:
        if configOptions[variantPopulation].user:
            configOptions[variantPopulation].value = set(configOptions[variantPopulation].value)
    seenVariants = set()
    for variantPopulation in ["pathological", "neutral", "predicting"]:
        if configOptions[variantPopulation].user:
            toDelete = seenVariants & configOptions[variantPopulation].value
            if toDelete:
                configOptions[variantPopulation].value = (
                    configOptions[variantPopulation].value - toDelete
                )
            seenVariants = seenVariants | configOptions[variantPopulation].value
            configOptions[variantPopulation].value = list(configOptions[variantPopulation].value)
    return configOptions


def retrieveVariants(configOptions):
    """
    Input: config options
    Function: retrieve variants from user, humsavar, sysMut.
              No repeated variants within a group and between groups. Patho neu pred priority.
              msa(DB/Build) variants will be retrieved later.
    Output: config options with values filled. Discarded if any variant is in bad format
    """

    # Initialize values
    for variantPopulation in ["pathological", "neutral", "predicting"]:
        if configOptions[variantPopulation].user:
            configOptions[variantPopulation].value = []

    # Retrieve User Variants
    for variantPopulation in ["pathological", "neutral", "predicting"]:
        if (
            configOptions[variantPopulation].user
            and variantPopulation + "DB" in configOptions[variantPopulation].user
        ):
            configOptions[variantPopulation].value.extend(
                retrieveUserVariants(
                    configOptions[variantPopulation + "DB"].user,
                    configOptions["protSeq"].user,
                )
            )

    # Retrieve Humsavar Variants
    humsavarPopulations = [
        variantPopulation
        for variantPopulation in ["pathological", "neutral"]
        if configOptions[variantPopulation].user
        and "humsavarDB" in configOptions[variantPopulation].user
    ]
    if humsavarPopulations:
        pathoHumsavar, neuHumsavar = retrieveHumsavarVariants(
            configOptions["humsavarDB"].user,
            configOptions["uniprotID"].user,
            configOptions["protSeq"].user,
        )
        if "pathological" in humsavarPopulations:
            configOptions["pathological"].value.extend(pathoHumsavar)
        if "neutral" in humsavarPopulations:
            configOptions["neutral"].value.extend(neuHumsavar)

    # Retrieve sysMut Variants
    if configOptions["predicting"].user and "sysMut" in configOptions["predicting"].user:
        configOptions["predicting"].value = retrieveSysMutVariants(configOptions["protSeq"].user)

    configOptions = pathoNeuPredPriority(configOptions)

    return configOptions


def retrieveAlignedQuery(msaDB, uniprotID, protSeq):
    """
    Input: msaDB, uniprotID, protSeqs
    Function: retrieve alignedQuery from msaDB and check that sequence is the same as user / blastDB
    Output: alignedQuery / sys.exit() if missing query protein in blast or is not equal as user / blastDB
    """

    for record in msaDB:
        if "UniRef100_" + uniprotID in record.id.split("|") or uniprotID in record.id.split("|"):
            alignedQuery = record
            if str(alignedQuery.seq).replace("-", "") != protSeq:
                sys.exit(
                    "ERROR: protein sequence from msaDB is not the same as protein sequence of user/blastProtDB"
                )
            return alignedQuery
    sys.exit("ERROR: missing " + uniprotID + " protein in msaDB")


def retrieveMSAFile(filename, uniprotID, protSeq):
    """
    Input: msa file, uniprotID, protSeq
    Function: retrieve msa and aligned query protein from file
    Output: msa / sys.exit() if missing query protein or sequence is not the same as blastProtDB/user
    """

    msaDB = list(Bio.SeqIO.parse(open(filename), "fasta"))
    if not msaDB:
        sys.exit("ERROR: msaDB empty")
    msaAa = "^(" + "|".join(list(extended_protein_letters + "-")) + ")+$"
    for rec in msaDB:
        if not re.search(msaAa, str(rec.seq)):
            sys.exit("ERROR: sequence of " + rec.id + " in msaDB contains a strange character")
    alignedQuery = retrieveAlignedQuery(msaDB, uniprotID, protSeq)
    return msaDB, alignedQuery


def retrieveMSANeutralVariants(configOptions):
    """
    Input: configOptions
    Function: retrieve neutral variants from MSA and apply patho>neu>pred priority
    Output: configOptions
    """

    alignedQuery = configOptions["alignedQuery"].value

    configOptions["neutralMSAID"].value = []
    alignedQueryLength = len(str(alignedQuery.seq).replace("-", ""))
    for index, msaSeq in enumerate(configOptions["msaDB"].value):
        if msaSeq.id != alignedQuery.id:
            identity = len(
                [
                    aa
                    for position, aa in enumerate(msaSeq.seq)
                    if aa == alignedQuery.seq[position] and alignedQuery.seq[position] != "-"
                ]
            ) / (float(len(str(msaSeq.seq).replace("-", "")) + alignedQueryLength) / 2)
            if configOptions["idmin"].value <= identity <= configOptions["idmax"].value:
                residue = 0
                for position, nAA in enumerate(alignedQuery.seq):
                    if nAA != "-":
                        residue += 1
                        if (
                            msaSeq.seq[position] != "-"
                            and msaSeq.seq[position] != nAA
                            and nAA + str(residue) + msaSeq.seq[position]
                            and checkNaturalAa(nAA)
                            and checkNaturalAa(msaSeq.seq[position])
                        ):
                            variant = nAA + str(residue) + msaSeq.seq[position]
                            configOptions["neutralMSAID"].value.append((variant, msaSeq.id))

    if configOptions["neutral"].user:
        configOptions["neutral"].value.extend([i for i, j in configOptions["neutralMSAID"].value])
        configOptions = pathoNeuPredPriority(configOptions)

    return configOptions


def loadConfigOptions():
    """
    Function: set up a list of available configuration options
    Output: list of option objects
    """

    configOptions = {
        # Input
        "uniprotID": option(
            "regex",
            options="^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}|[OPQ][0-9][A-Z0-9]{3}[0-9]-[0-9]{1,2}|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}-[0-9]{1,2})$",
        ),
        "protSeq": option(),
        # Function: MSA only
        "msaDB": option("path"),
        "msaBuild": option(
            "boolean",
            dependencies="blastPath_blastProtDB_numIter_evalue_minSeqIden_muscleAlg",
        ),
        # Function: ARFF too
        "pathological": option("list", options="humsavarDB,pathologicalDB", dependencies="user"),
        "neutral": option(
            "list", options="humsavarDB,neutralDB,msaDB,msaBuild", dependencies="user"
        ),
        "predicting": option("list", options="predictingDB,sysMut", dependencies="user"),
        # Output
        "output": option("path"),
        "write": option("list", options="arff,msa,neutresMsaID"),
        # MSA Build Options
        "minSeqIden": option("float", options="0-100"),
        "numIter": option("int", options="2-2"),
        "evalue": option("float", options="0-1"),
        # Model neutral variants from MSA
        "idmin": option("float", options="0-1"),
        "idmax": option("float", options="0-1"),
        # ARFF Features
        "impRes": option("boolean", dependencies="uniprotImpResDB|moreImpResDB"),
        "vdwVolume": option("boolean", dependencies="vdwVolumeDB"),
        "hydrophobicity": option("boolean", dependencies="hydrophobicityDB"),
        "substitutionMatrix": option("boolean", dependencies="substitutionMatrixDB"),
        "pssm-native": option("boolean"),
        "pssm-mutated": option("boolean"),
        "entropy": option("boolean"),
        "freqGaps": option("boolean"),
        "entropyAvgWindow": option("boolean", dependencies="entropy"),
        "relEntropy": option("boolean"),
        "fractionSimilarity": option("boolean"),
        # Databases
        "blastProtDB": option("path"),
        "humsavarDB": option("path"),
        "pathologicalDB": option("path"),
        "neutralDB": option("path"),
        "predictingDB": option("path"),
        "uniprotImpResDB": option("path"),
        "moreImpResDB": option("path"),
        "vdwVolumeDB": option("path"),
        "hydrophobicityDB": option("path"),
        "substitutionMatrixDB": option(
            "regex", options="^(" + "|".join(substitution_matrices.load()) + ")$"
        ),
        # Softwares
        "blastPath": option("path"),
        "muscleAlg": option("path"),
        # More
        "absoluteFlag": option("boolean"),
        "alignedQuery": option(),
        "neutralMSAID": option(),
        "variants": option(),
        "tag": option(),
        "impResDB": option(),
        "msaDBT": option(),
        "config": option(),
    }

    return configOptions


def parseConfig(configFile, configOptions):
    """
    Input: configuration file, configurationOptions
    Function: parse user's configuration from a file and add them to configurationOptions
    Output: configurationOptions
    """

    with open(configFile) as fh:
        for lineNum, line in enumerate(fh.readlines()):
            line = line.rstrip()
            if line and not line.startswith("#"):
                if not "=" in line:
                    sys.exit(
                        "ERROR: line "
                        + str(lineNum)
                        + " of configuration file "
                        + configFile
                        + ' has a bad format. May be "=" is missing?'
                    )
                line = line.replace(" ", "").split("=")
                line = [line[0], "=".join(line[1:])]
                if line[0] not in configOptions.keys():
                    sys.exit(
                        "ERROR: line "
                        + str(lineNum)
                        + " of configuration file "
                        + configFile
                        + " has a nonexistent option: "
                        + line[0]
                        + ". Available options are: \n"
                        + "\n".join(sorted(configOptions.keys()))
                    )
                if line[1] != "False":
                    configOptions[line[0]].user = line[1]
    configOptions["config"].user = configFile
    return configOptions


def checkMandatory(configOptions):
    """
    Input: configurationOptions
    Function: check if mandatory options were introduced by the user
    Output: sys.exit() if any mandatory/dependency option is missing
    """

    # Check Mandatory: first level AND, second level OR
    mandatory = "uniprotID_pathological|neutral|predicting|msaBuild_msaDB|msaBuild_write_output"
    for options in mandatory.split("_"):
        if not any([configOptions[option].user for option in options.split("|")]):
            sys.exit("ERROR: missing mandatory options: " + options.replace("|", " or "))


def checkDependencies(configOptions):
    """
    Input: configurationOptions
    Function: check if dependencies derived from the user selected options are satisfied
    Output: sys.exit() if any mandatory/dependency option is missing
    """

    # Check Dependencies: first level OR, second level AND
    def checkDependenciesOption(deps, configOptions):
        for dependencies in deps.split("|"):
            if all([configOptions[dependency].user for dependency in dependencies.split("_")]):
                return True
        sys.exit(
            "ERROR: missing dependency options: " + deps.replace("_", " and ").replace("|", " or ")
        )

    # Check special dependencies
    if not configOptions["protSeq"].user:
        checkDependenciesOption("blastProtDB_blastPath", configOptions)
    if "arff" in configOptions["write"].user:
        checkDependenciesOption("pathological|neutral|predicting", configOptions)
    if "neutresMsaID" in configOptions["write"].user or (
        configOptions["neutral"].user
        and (
            "msaDB" in configOptions["neutral"].user or "msaBuild" in configOptions["neutral"].user
        )
    ):
        checkDependenciesOption("idmin_idmax", configOptions)
        if configOptions["idmin"].user > configOptions["idmax"].user:
            sys.exit("ERROR: idmin should be smaller or equal than idmax")
    if (
        "neutresMsaID" in configOptions["write"].user
        and "arff" in configOptions["write"].user
        and not (
            configOptions["neutral"].user
            and (
                "msaDB" in configOptions["neutral"].user
                or "msaBuild" in configOptions["neutral"].user
            )
        )
    ):
        sys.exit("ERROR: to calculate neutresMsaID is necessary msaDB or msaBuild in neutral.")

    # Check config dependencies
    dependencyExceptions = ["sysMut"]
    for option in configOptions.values():
        if option.user and option.dependencies:
            if option.dependencies == "user":
                for dependency in option.user:
                    if dependency not in dependencyExceptions:
                        checkDependenciesOption(dependency, configOptions)
            else:
                checkDependenciesOption(option.dependencies, configOptions)


def validateOptions(configOptions):
    """
    Input: configurationOptions
    Function: validate user options and compute what is needed: retrieve data from files, search protein sequence
    Output: configurationOptions
    """

    # Validate user options
    for name, option in configOptions.items():
        if option.user:
            if option.type == "boolean":
                if option.user != "True":
                    sys.exit(f"ERROR: {name} option must be True or False")
                else:
                    configOptions[name].user = True
            elif option.type == "path":
                if name == "output":
                    if os.path.isdir(option.user):
                        if not option.user.endswith("/"):
                            option.user += "/"
                        option.user += configOptions["uniprotID"].user
                        continue
                    elif os.path.isdir(os.path.dirname(option.user)):
                        continue
                    else:
                        sys.exit(f"ERROR: {option.user} path does not exists")
                if not os.path.exists(option.user):
                    if name == "blastProtDB":
                        if os.path.exists(f"{option.user}.00.phr"):
                            continue
                    sys.exit(f"ERROR: {option.user} path does not exists")
            elif option.type == "float":
                try:
                    configOptions[name].user = float(option.user)
                except ValueError:
                    sys.exit(f"ERROR: {option.user} is not a number")
                minValue, maxValue = option.options.split("-")
                if not float(minValue) <= option.user <= float(maxValue):
                    sys.exit(f"ERROR: {option.user} should be between {minValue}-{maxValue}")
            elif option.type == "int":
                try:
                    configOptions[name].user = int(option.user)
                except ValueError:
                    sys.exit(f"ERROR: {option.user} is not a number")
                minValue, maxValue = option.options.split("-")
                if not int(minValue) <= option.user <= int(maxValue):
                    sys.exit(f"ERROR: {option.user} should be between {minValue}-{maxValue}")
            elif option.type == "regex":
                if not re.search(option.options, option.user):
                    sys.exit(f"ERROR: {name} options with value {option.user} is misspelled")
            elif option.type == "list":
                for item in option.user.split(","):
                    if item not in option.options.split(","):
                        sys.exit(f"ERROR: {item} option not available")
                configOptions[name].user = option.user.split(",")

    return configOptions


def computeOptions(configOptions):
    """
    Input: configurationOptions
    Function: validate user options and compute what is needed: retrieve data from files, search protein sequence
    Output: configurationOptions
    """

    # Compute values
    # Retrieve Mandatory
    if not configOptions["protSeq"].user:
        configOptions["protSeq"].user = retrieveProtSeq(
            configOptions["uniprotID"].user,
            configOptions["blastProtDB"].user,
            configOptions["blastPath"].user,
        )
    if not checkNaturalAa(configOptions["protSeq"].user):
        sys.exit("ERROR: protein sequence has an unnatural amino acid")

    # Retrieve optionals
    if configOptions["substitutionMatrixDB"].user:
        configOptions["substitutionMatrixDB"].value = {}
        for oriaa in "ACDEFGHIKLMNPQRSTVWY":
            for mutaa in "ACDEFGHIKLMNPQRSTVWY":
                if oriaa != mutaa:
                    configOptions["substitutionMatrixDB"].value[
                        (oriaa, mutaa)
                    ] = substitution_matrices.load(configOptions["substitutionMatrixDB"].user)[
                        oriaa
                    ][
                        mutaa
                    ]
    if configOptions["vdwVolumeDB"].user:
        configOptions["vdwVolumeDB"].value = file2dict(configOptions["vdwVolumeDB"].user)
    if configOptions["hydrophobicityDB"].user:
        configOptions["hydrophobicityDB"].value = file2dict(configOptions["hydrophobicityDB"].user)
    if configOptions["uniprotImpResDB"].user or configOptions["moreImpResDB"].user:
        configOptions["impResDB"].value = []
        if configOptions["uniprotImpResDB"].user:
            configOptions["impResDB"].value.extend(
                retrieveImpRes(
                    configOptions["uniprotImpResDB"].user,
                    configOptions["uniprotID"].user,
                    len(configOptions["protSeq"].user),
                )
            )
        if configOptions["moreImpResDB"].user:
            configOptions["impResDB"].value.extend(
                retrieveImpRes(
                    configOptions["moreImpResDB"].user,
                    "",
                    len(configOptions["protSeq"].user),
                )
            )
    if "arff" in configOptions["write"].user or "neutresMsaID" in configOptions["write"].user:
        configOptions = retrieveVariants(configOptions)
    if configOptions["msaDB"].user:
        (configOptions["msaDB"].value, configOptions["alignedQuery"].value,) = retrieveMSAFile(
            configOptions["msaDB"].user,
            configOptions["uniprotID"].user,
            configOptions["protSeq"].user,
        )

    # Copy user to values options if need it
    for option in configOptions:
        if configOptions[option].user and not configOptions[option].value:
            if option not in ["neutral", "pathological", "predicting"]:
                configOptions[option].value = configOptions[option].user

    return configOptions


def version():
    """
    Function: Obtain git version of the program
    Output: program version
    """

    gitPath = f"{os.path.dirname(os.path.abspath(__file__))}/.git"
    gitVersion = shell(f"git --git-dir={gitPath} rev-list --count HEAD").rstrip("\n")
    return gitVersion


def commandLine():
    """
    Input: command line
    Function: parse configuration file or show usage, help or version of software
    Output: sys.exit() if help need it
    """

    def help():
        return f"usage: python PatMut.py <config file>\noptions\n-v, -version\tPatMut v{version()}\n-h, --help"

    # Parse options
    if len(sys.argv) != 2 or sys.argv[1] in ["-h", "--help"]:
        sys.exit(help())
    elif sys.argv[1] in ["-v", "--version"]:
        sys.exit("PatMut v" + version())
    # Retrieve config file
    elif os.path.isfile(sys.argv[1]):
        return sys.argv[1]
    elif not os.path.exists(sys.argv[1]):
        sys.exit(sys.argv[1] + " path does not exists")
    # What could else happen?
    else:
        sys.exit(help())


def configuration(configFile):
    """
    Input: configuration file
    Function: set up configuration options
    Output: configuration options
    """

    # Load patmut options
    configOptions = loadConfigOptions()

    # Parse user options
    configOptions = parseConfig(configFile, configOptions)

    # Check mandatory options & option's dependencies
    checkMandatory(configOptions)

    configOptions = validateOptions(configOptions)

    checkDependencies(configOptions)

    # Validate user's optionValues and compute final optionValues
    configOptions = computeOptions(configOptions)

    return configOptions


def config():
    """
    Input: command line
    Function: main
    Output: configOptions
    """

    configFile = commandLine()
    configOptions = configuration(configFile)

    return configOptions
