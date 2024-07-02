# PatMut
PatMut extracts physicochemical and sequence features of missense variants.

# Requirements

You need to download the following datasets:

- **UniRef100 fasta**: UniProt database of clustered sets of protein sequences. You need to download it from the [UniProt FTP site](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/) and convert it to a BLAST database with the command `makeblastdb`.
- **Humsavar**: UniProt database of curated pathogenic and benign variants from the literature. It is already provided in the database directory. You can download it from the [UniProt FTP site](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/humsavar.txt).
- **Important residues**: database of important residues. The one provided in the database directory is based on the information in UniProt. You can create your own following the format of one protein per row with its UniProt identifier, a tab, and the list of important residues separeted by a comma. 
- **Amino acid volumes**: database of amino acid's volumes. The one provided in the database directory is based on theVan der Waals amino acid volumes of Bondi A. Van der Waals volumes and radii. J Phys Chem 1964;68:441–451. You can create your own following the format of one amino acid per row with the amino acid identifier, a tab, and the volume.
- **Amino acid hydrophobicities**: database of amino acid's hydrophobicities. The one provided in the database directory is based on the water/octanol free energy measurements of Fauchere JL, Pliska V. Hydrophobic parameters of amino acid side chains from the partitioning of N-acetyl-amino acid amides. Eur J Med Chem 1983;18:369–375. You can create your own following the format of one amino acid per row with the amino acid identifier, a tab, and the hydrophobicity.

Additionally, you need to install the following software:

- **BLAST+ executables**: download and install BLAST+ executables from [NCBI BLAST FTP site](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
- **MUSCLE**: download and install MUSCLE from [MUSCLE website](https://www.drive5.com/muscle/).

# Execute PatMut

## Create the configuration file

The configuration file is a set of key-value pairs of parameter=value. There are the following parameters:

- **QUERY**
  - `uniprotID`: the UniProt identifier of the protein of interest.
  - `blastPath`: path to the directory where the BLAST executable is located.
  - `blastProtDB`: path to the UniRef100 BLAST protein database.  

- **RESULTS**
  - `output`: path to the directory where all the results are saved. 
  - `write`: list of results that will be saved separated by a comma. The possible results are:
    - `arff`: list of features of the variants of the protein in [arff format](https://waikato.github.io/weka-wiki/formats_and_processing/arff_stable/#). The features are defined in the PROPERTIES section and the variants are configured in the MUTATION SET section.
    - `msa`: MSA of the protein of interest and other related to it in FASTA format. 
    - `neutresMsaID`: list of neutral variants identified in the MSA.


- **MUTATION SET**
  - `pathological`: source of pathological variants. These variants have a tag of 1. Possible values are `humsavarDB` or `pathologicalDB`.
  - `neutral`: source of neutral variants. These variants have a tag of 0. Possible values are `humsavarDB`, `neutralDB`, `msaDB` or `msaBuild`.
  - `predicting`: source of uncategorized variants. These variants have a tag of ?. Possible values are `predictingDB` or `sysMut`.

- **PROPERTIES**

    List of physicochemical and sequence features that needs to be calculated for the variants. The key is the name of the feature and the value is the keyword `True`. Possible features are `vdwVolume`, `hydrophobicity`, `substitutionMatrix`, `pssm-native`, `entropy`, `impRes`. Additoinally, there is the parameter `absoluteFlag` to indicate if the hydrophobicity and volume change should be an absolute value or not. Possible values are `True` or `False`.

- **DATABASES**
  
    Paths to the database files of `humsavarDB`, `uniprotImpResDB`, `vdwVolumeDB`, `hydrophobicityDB`. Additionally, for the `substitutionMatrixDB` it needs to be indicated the name of the substitution matrix, such as `BLOSUM62`.

- **BUILD MSA**
  
    Section only specified when we want to build the MSA of a protein of interest. 
    
    Building the MSA is a slow process, so we advise to save the generated MSA in case we want to run again this protein. To use an already generated MSA we have to use the `msaDB` parameter and specify the path to the MSA.
    
    Parameters of Build MSA section are:
    -  `msaBuild`: indicates if we want to build the MSA. It needs to be set to `True`.
    - `numIter`: number of iterations of the BLAST search.
    - `evalue`: minimal e-value of the BLAST search.
    - `minSeqIden`: minimal sequence identity of the BLAST search.
    - `muscleAlg`: path to directory where the MUSCLE alignment tool is located.

- **Modelling Neutral Variants from MSA**
    
    When the source of neutral variants is msaDB or msaBuild, the Modelling Neutral Variants from MSA is executed. This section has two parameters `idmin` and `idmax` that accept a numerical value from 0 - 1. These parameters define the range of sequence identities accepted to model the neutral variants from the provided / built MSA.
 

## Run PatMut

To run PatMut execute in the terminal the python script PatMut.py with the desired configuration file, for example:

```
python bin/PatMut.py demo/example_provide_msa.config
```