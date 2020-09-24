<<<<<<< HEAD
# LaanLab-SATAY-DataAnalysis

This repository contains codes, data and workflows for data analysis regarding SATAY experiments.

## [Projects tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/projects/1) and [Issues tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/issues)

For comments, ToDo lists and general progress updates on the projects, please notify this in the [Projects tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/projects/1).
Issues and questions about codes can be published in the [Issues tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/issues).

## [Docs](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/docs)

Files with detailed notes about analysis, datafiles, software and workflows (including an installation guide for the Linux workflow for SATAY analysis).

## [Data files](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Data_Files)

Currently the following data files and folder are present:

1. [Cerevisiae_EssentialGenes_list_1.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cervisiae_EssentialGenes_List_1.txt): This is a list of known essential genes with systematic naming format.
2. [Cerevisiae_EssentialGenes_List_2.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cervisiae_EssentialGenes_List_2.txt): This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes.
3. [Cerevisiae_AllEssentialGenes_list.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_AllEssentialGenes_List.txt): This is a combination of the previous two files where all the duplicates have been removed and the standard yeast gene naming convention has been applied.
4. [Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Yeast_Protein_Names.txt): This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
5. [S288C_reference_sequence_R64-2-1_20150113.fsa](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S288C_reference_sequence_R64-2-1_20150113.fsa): Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
6. [Saccharomyces_cerevisiae.R64-1-1.99.gff3](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Saccharomyces_cerevisiae.R64-1-1.99.gff3): Information file for all CDS in the *S.Cerevisiae* genome in .gff format.
7. [S_Cerevisiae_protein_designation_name_full_genome.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S_Cerevisiae_protein_designation_name_full_genome.txt): Includes all gene names present in *S.Cerevisiae* in designation naming convention.
8. [S_Cerevisiae_protein_oln_name_full_genome.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S_Cerevisiae_protein_oln_name_full_genome.txt): Includes all gene names present in *S.Cerevisiae* in oln naming convention.
9. [SGD_features.tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/SGD_features.tab): Includes all features (e.g. genes, ars, telomeres etc.) present in the *S.Cerevisiae* genome (this file has a [readme](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/SGD_features.README)).

### [test_data](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Data_Files/test_data)

This folder (located in the Data files folder) contains six files that are the result of the processing by the python code for transposon analysis (see the folder 'python_TransposonMapping; transposonsmapping_satay.py' (works only in Linux, for Windows try Matlab_TransposonMapping; tn_and_reads_per_gene.m, which generates only the first three files in the below list)).

1. [... .wig](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam.wig): Contains the number of reads per insertion location in wig format.
2. [... .bed](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam.bed): Contains the number of reads oer insertion location in bed format. This file is similar to the .wig file, but in the .bed file the number of reads are given by (reads*20)+100
3. [... _pergene.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam_pergene.txt): Contain the number of insertions and number of reads per gene. For this an overview is created of most known genes (6600 genes are considered).
4. [... _peressential.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam_peressential.txt): Contain the number of insertions and number of reads per essential gene. For this an overview is created for the annotated essential genes (1186 genes are considered).
5. [... _pergene_insertions.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam_pergene_insertions.txt): Contain the exact locations and reads for each insertion in each gene (6600 genes are considered).
6. [... _peressential_insertions.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/ERR1533148_trimmed.sorted.bam_peressential_insertions.txt): Contain the exact locations and reads for each insertion in each annotated essential gene (6600 genes are considered).

## [Python scripts](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/python_scripts) and [python modules](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/python_modules)

Data analysis created with python are located in the folder for python scripts.
Python functions are placed the python modules folder.
Each python code should include information on how and when to use it.
Some codes also come a python notebook format using the .ipynb extension.

## [Matlab_TransposonMapping](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Matlab_TransposonMapping)

Currently the following Matlab codes and data files are present:

1. [tn_and_reads_per_gene.m](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/tn_and_reads_per_gene.m): This inputs a .bam file and outputs the number of transposons and reads per gene.
2. [names.mat](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/names.mat): Matlab data file with all genes names. This file is required for running the matlab code 1.
3. [yeastGFF.mat](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/yeastGFF.mat): Matlab data file with information about the genes. This file is required for running [the matlab code](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/tn_and_reads_per_gene.m).

## [Python_TransposonMapping](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Python_TransposonMapping)

This python is used in the Linux workflow for analyzing SATAY data.
The file [transposonmapping_satay.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Python_TransposonMapping/transposonmapping_satay.py) includes the actual code and the folder python_modules include some modules necessary to run it.
This folder also include [processing_workflow shell script](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Python_TransposonMapping/processing_workflow.sh) that is used to call all software tools (including this python script) within Linux.
For setting up a virtual environment including all software tools required for running this workflow, see the [installation guide](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/docs/Installation_Guide_SATAY_Analysis_Software.pdf) in the [docs](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/docs) folder on this Github page.
Note: The python code only works in Linux operating systems. For running a similar code on Windows, try the [tn_and_reads_per_gene.m matlab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/tn_and_reads_per_gene.m) code, which generates a subset of the files that are generated in the python code (see also the 'test_data' section in this readme file).

## License
[LaanLab for Bionanoscience, Delft university of Technology](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/bionanoscience/research/research-labs/liedewij-laan-lab/research-projects/evolvability-and-modularity-of-essential-functions-in-budding-yeast/)

*Last updated: September 7, 2020*
=======
# LaanLab-SATAY-DataAnalysis (dev_gregory)

This repository contains codes, data and workflows for data analysis regarding SATAY experiments, created by Gregory van Beek.

For an example test dataset, see the Data Files folder in the master-branch.

## docs

Contains files for discussing analysis steps for SATAY data and installation guides for the required software.
For detailed information about all the software tools used and the experimental approaches for satay experiments, see the 'satay_analysis_notes'.

## Matlab_TransposonMapping

Contains files for transposon mapping based using .bam files.
This code is originally created by [the Kornmann-lab](https://sites.google.com/site/satayusers/complete-protocol/bioinformatics-analysis/matlab-script).
Currently the following Matlab codes and data files are present:

1. tn_and_reads_per_gene.m: This inputs a .bam file and outputs the number of transposons and reads per gene.
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code 1.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code 1.

## Python_TransposonMapping

'UNDER DEVELOPMENT'.
Contains files for transposon mapping using .bam files.
These python codes are based on the Matlab codes developed by the Kornmann lab (see codes in Matlab_TransposonMapping).
The python codes run only in Linux and require the package [pysam](https://pysam.readthedocs.io/en/latest/index.html) installed from Bioconda.

## Python_scripts

Contains python codes for different kind of analysis. The most important codes are also given as python notebooks with more detailed explanation and interpretation of the results.
For each python script a .ipynb file is created that includes the same code as the .py file with the same name.
These .ipynb files contain more explaination about the codes and the interpretation of the results.

### Python_modules

This folder contains codes for extracting information about genes from online data files. All required text files are present in the Data_Files folder on this branch.
The following files are present:

1. chromosome_and_gene_positions.py: This contains the functions `chromosome_positions` and `gene_positions` which outputs information about the location and lengths of the chromosomes and genes. Input is the gff file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'. Also contains the function `chromosome_roman_to_arabic` which translates the roman numerals to arabic numerals or vice versa for the 16 chromosomes.
2. chromosome_names_in_files.py: Contains the functions `chromosome_name_bedfile` and `chromosome_name_wigfile` which can be used for getting the chromosome names used in the bed file or wig file, respectively. This is needed since the chromosome names in the files are not always given as roman numerals, but this is often expected in other codes.
3. essential_genes_names.py: This contains the function `list_known_essentials`. It combines and outputs the essential genes from the files 'Cervisiae_EssentialGenes_List_1.txt' and 'Cervisiae_EssentialGenes_List_2.txt' in one list.
4. gene_length.py: This contains the functions `gene_length_bp` and `gene_length_aa`, both which requires the file 'Yeast_Protein_Names.txt'. It extracts the gene length in terms of basepairs and amino-acids, respectively.
5. gene_names.py: This contains the functions `list_gene_names` and `gene_aliases` and both require the file 'Yeast_Protein_Names.txt' as input. The function 'list_gene_names' creates a list of all gene names, including aliases and different naming conventions. The function 'gene_aliases' creates dictionaries with keys the gene names according to the systematic naming convention (oln) and values are the names according to the standard naming conventions (including all potential aliases of the genes), the gene_id for SGD and the gene_id for swiss-prot. Input is the gff file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'.
6. statistics_perchromosome.py: Contains the function `chromosome_insertion_periodicity` which prints some statistical values for transposon insertion per chromosome. Output is the distances between insertions given in terms of basepairs.

### Data_Files

Currently the following data files are present:

1. Cerevisiae_EssentialGenes_list_1.txt: This is a list of known essential genes with systematic naming format.
2. Cerevisiae_EssentialGenes_List_2.txt: This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes.
3. Yeast_Protein_Names.txt: This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
4. S288C_reference_sequence_R64-2-1_20150113.fsa: Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
5. Saccharomyces_cerevisiae.R64-1-1.99.gff3: Contains an overview of all the Coding DNA sequences in the yeast genome including information about these sequences.

*Last updated: June 09, 2020*
>>>>>>> 02723fc28becaa0428d64ab6978185f5da1aa080
