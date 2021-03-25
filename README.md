# LaanLab-SATAY-DataAnalysis (dev_gregory)
test
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
2. names.mat: Matlab data file with all genes names. This file is required for running the matlab code one.
3. yeastGFF.mat: Matlab data file with information about the genes. This file is required for running the matlab code one.

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
