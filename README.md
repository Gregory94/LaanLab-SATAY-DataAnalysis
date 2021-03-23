
# LaanLab-SATAY-DataAnalysis

This repository contains codes, data and workflows for data analysis regarding SATAY experiments.
This branch is the old master branch of this repository.

## License

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
This work is licensed under Apache 2.0 . 
The 2.0 version of the Apache License, approved by the ASF in 2004, helps us achieve our goal of providing reliable and long-lived software products through collaborative open source software development.

## [Projects tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/projects/1) and [Issues tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/issues)

For comments, ToDo lists and general progress updates on the projects, please notify this in the [Projects tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/projects/1).
Issues and questions about codes can be published in the [Issues tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/issues).

## [Docs](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/docs)

Files with detailed notes about analysis, datafiles, software and workflows (including an installation guide for the Linux workflow for SATAY analysis).

## [Data files](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Data_Files)

Currently the following data files and folder are present:

1. [Cerevisiae_EssentialGenes_list_1.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_EssentialGenes_List_1.txt): This is a list of known essential genes with systematic naming format.
2. [Cerevisiae_EssentialGenes_List_2.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_EssentialGenes_List_2.txt): This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes.
3. [Cerevisiae_AllEssentialGenes_list.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_AllEssentialGenes_List.txt): This is a combination of the previous two files where all the duplicates have been removed and the standard yeast gene naming convention has been applied.
4. [Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Yeast_Protein_Names.txt): This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
5. [S288C_reference_sequence_R64-2-1_20150113.fsa](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S288C_reference_sequence_R64-2-1_20150113.fsa): Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
6. [Saccharomyces_cerevisiae.R64-1-1.99.gff3](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Saccharomyces_cerevisiae.R64-1-1.99.gff3): Information file for all CDS in the *S.Cerevisiae* genome in .gff format.
7. [S_Cerevisiae_protein_designation_name_full_genome.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S_Cerevisiae_protein_designation_name_full_genome.txt): Includes all gene names present in *S.Cerevisiae* in designation naming convention.
8. [S_Cerevisiae_protein_oln_name_full_genome.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/S_Cerevisiae_protein_oln_name_full_genome.txt): Includes all gene names present in *S.Cerevisiae* in oln naming convention.
9. [SGD_features.tab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/SGD_features.tab): Includes all features (e.g. genes, ars, telomeres etc.) present in the *S.Cerevisiae* genome (this file has a [readme](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/SGD_features.README)).

### [test_data](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Data_Files/test_data)

This folder (located in the Data files folder) contains four subfolders containing each seven files that are the result of the processing by the python code for transposon analysis (see the folder '[python_TransposonMapping](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Python_TransposonMapping); [transposonsmapping_satay.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Python_TransposonMapping/transposonmapping_satay.py)' (works only in Linux, for Windows try Matlab_TransposonMapping; tn_and_reads_per_gene.m, which generates only the first three files in the below list)).
Each subfolder contains information for a specific dataset.
The location of the raw data can be found in the link stored in [KornmannLab_dataset_Link.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/KornmannLab_Dataset_Link.txt).
The links in the following list direct to the dataset in the folder KornmannLab_dataset_WT1, but similar files are found in all dataset folders.

1. [... .wig](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam.wig): Contains the number of reads per insertion location in wiggle format.
2. [... .bed](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam.bed): Contains the number of reads per insertion location in Browser Extensible Data format. This file contains the same information as the .wig file, but in the .bed file the number of reads are given by (reads * 20) + 100. In the bed-format, the insertions at the same location but with a different orientation are stored as individual insertions, opposed to the wig-files where insertions at the same location with different orientations are added up.
3. [... _pergene.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam_pergene.txt): Contain the number of insertions and number of reads per gene. For this an overview is created of most known genes (6600 genes are considered).
4. [... _peressential.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam_peressential.txt): Contain the number of insertions and number of reads per essential gene. For this an overview is created for the annotated essential genes (1186 genes are considered).
5. [... _pergene_insertions.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam_pergene_insertions.txt): Contain the exact locations and reads for each insertion in each gene (6600 genes are considered).
6. [... _peressential_insertions.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_trimmed.sorted.bam_peressential_insertions.txt): Contain the exact locations and reads for each insertion in each annotated essential gene (6600 genes are considered).
7. [... _log.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/test_data/Kornmannlab_dataset_WT1/ERR1533147_log.txt): Log file containing the details for the processing procedure as completed with the workflow as described in the [Python_TransposonMapping folder](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Python_TransposonMapping).

## [Python scripts](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/Python_scripts) and [python modules](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/python_modules)

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
This folder also include [processing_workflow shell script](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/satay_processing_workflow.sh) that is used to call all software tools (including this python script) within Linux.
For setting up a virtual environment including all software tools required for running this workflow, see the [installation guide](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/docs/Installation_Guide_SATAY_Analysis_Software.pdf) in the [docs](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/docs) folder on this Github page.
Note: The python code only works in Linux operating systems. For running a similar code on Windows, try the [tn_and_reads_per_gene.m matlab](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Matlab_TransposonMapping/tn_and_reads_per_gene.m) code, which generates a subset of the files that are generated in the python code (see also the 'test_data' section in this readme file).

## License
[LaanLab for Bionanoscience, Delft university of Technology](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/bionanoscience/research/research-labs/liedewij-laan-lab/research-projects/evolvability-and-modularity-of-essential-functions-in-budding-yeast/)

*Last updated: January 7, 2021*
