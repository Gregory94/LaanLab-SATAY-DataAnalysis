# LaanLab-SATAY-DataAnalysis
This repository contains codes, data and workflows for data analysis regarding SATAY experiments, created by Gregory van Beek.

## Data files
Currently the following data files are present:

1. Cerevisiae_EssentialGenes_list_1.txt: This is a list of known essential genes with systematic naming format.
2. Cerevisiae_EssentialGenes_List_2.txt: This is a list of known essential genes with common naming format. Some genes may occur only in one file, so it is recommended to use both files simultaneously to have a complete list of the known essential genes. 
3. Yeast_Protein_Names.txt: This is a list that includes all genes with both naming convention and their ID's. It also include the length of the corresponding proteins in terms of amino acids.
4. S288C_reference_sequence_R64-2-1_20150113.fsa: Reference sequence for wild type cells of *S.Cerevisiae* from the S288C strain.
5. E-MTAB-4885.WT1.bam_pergene.txt: This file is an example of how the result looks like after processing and can be used as a test file. It contains a list of all analyzed genes with the respective read and transposon counts. This file can be input in the statistics_pergene.py code.

## docs
Contains files for discussing analysis steps for SATAY data and installation guides for the required software.

## Matlab_TransposonMapping
Contains files for transposon mapping based on .bam files. This code is originally created by the Kornmann-lab.
Currently the following Matlab codes and data files are present:

1. tn_and_reads_per_gene.m: This inputs a .bam file and outputs the number of transposons and reads per gene. The output of this code can be used for python notebook 1.
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code 1.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code 1.

## Python_modules
This folder contains codes for extracting information about genes from online data files. All required text files are present in the docs folder on this branch.
The following files are present:

1. chromosome_and_gene_positions.py: This contains the functions `chromosome_positions` and `gene_positions` which outputs information about the location and lengths of the chromosomes and genes. Input is the gff file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'.
2. essential_genes_names.py: This contains the function `list_known_essentials`. It combines and outputs the essential genes from the files 'Cervisiae_EssentialGenes_List_1.txt' and 'Cervisiae_EssentialGenes_List_2.txt' in one list.
3. gene_length.py: This contains the functions `gene_length_bp` and `gene_length_aa`, both which requires the file 'Yeast_Protein_Names.txt'. It extracts the gene length in terms of basepairs and amino-acids, respectively.
4. gene_names.py: This contains the functions `list_gene_names` and `gene_aliases` and both require the file 'Yeast_Protein_Names.txt' as input. The function 'list_gene_names' creates a list of all gene names, including aliases and different naming conventions. The function 'gene_aliases' creates dictionaries with keys the gene names according to the systematic naming convention (oln) and values are the names according to the standard naming conventions (including all potential aliases of the genes), the gene_id for SGD and the gene_id for swiss-prot. Input is the gff file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3'.

## Other files
Other files are python files for coding that are not yet suitable for the master branch.

*Last updated: April 17, 2020*
