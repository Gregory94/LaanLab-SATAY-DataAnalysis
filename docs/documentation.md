# Documentation for processing SATAY <!-- omit in toc -->

- [Introduction](#introduction)
- [File types](#file-types)
  - [fastq](#fastq)
  - [sam & bam](#sam--bam)
  - [wig](#wig)
  - [bed](#bed)
  - [pergene.txt & peressential.txt](#pergenetxt--peressentialtxt)
  - [pergene_insertions.txt & peressential_insertions.txt](#pergene_insertionstxt--peressential_insertionstxt)
- [Software - Processing](#software---processing)
  - [satay.sh](#sataysh)
- [Software - analysis](#software---analysis)
  - [python scripts](#python-scripts)
    - [strip_redundant_insertions.py](#strip_redundant_insertionspy)
    - [genomicfeatures_dataframe.py](#genomicfeatures_dataframepy)
    - [scatterplot_genes.py](#scatterplot_genespy)
    - [transposonread_profileplot.py](#transposonread_profileplotpy)
    - [transposonread_profileplot_genome.py](#transposonread_profileplot_genomepy)
    - [volcanoplot.py](#volcanoplotpy)
    - [create_essentialgenes_list.py](#create_essentialgenes_listpy)
    - [split_wigfiles.py](#split_wigfilespy)
  - [python modules](#python-modules)
    - [chromosome_and_gene_positions.py](#chromosome_and_gene_positionspy)
    - [chromosome_names_in_files.py](#chromosome_names_in_filespy)
    - [dataframe_from_pergene.py](#dataframe_from_pergenepy)
    - [essential_genes_names.py](#essential_genes_namespy)
    - [gene_length.py](#gene_lengthpy)
    - [gene_names.py](#gene_namespy)
    - [gene_tn_insertions.py](#gene_tn_insertionspy)
    - [insertions_count.py](#insertions_countpy)
    - [mapped_reads.py](#mapped_readspy)
    - [read_sgdfeatures.py](#read_sgdfeaturespy)
    - [statistics_perchromosome.py](#statistics_perchromosomepy)
  - [Other tools](#other-tools)
    - [IGV](#igv)
    - [genome browser](#genome-browser)
- [Outlook](#outlook)
- [How to use the Linux desktop](#how-to-use-the-linux-desktop)
- [Appendices](#appendices)
  - [PHRED tables (base33)](#phred-tables-base33)
  - [PHRED table (base64)](#phred-table-base64)

This documentation gives a complete overview for the processing of the data from SAturated Transposon Analysis in Yeast (SATAY).
It includes a short introduction to SATAY and a detailed discussion on to perform the processing from the raw sequencing data to the postprocessing and checking of the results.

For the processing a pipeline is created using Bash and Python.
The workflow and the python codes can be found at [github.com/Gregory94/LaanLab-SATAY-DataAnalysis](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/satay_processing).
More information about satay analysis and experimental protocols can be found at the [satayusers website from the Kornmann lab](https://sites.google.com/site/satayusers/ "satayusers website") or, for more questions, visit the [satayusers forum](https://groups.google.com/g/satayusers "satayusers forum").

> Date last update: 10-03-2021
> 
> Contact: Gregory van Beek (G.M.vanBeek@tudelft.nl)
>  
> [Laanlab, Delft University of Technology](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/bionanoscience/research/research-labs/liedewij-laan-lab/research-projects/evolvability-and-modularity-of-essential-functions-in-budding-yeast "LaanLab TUDelft")

## Introduction

## File types

During the processing, different file types are being used and created.
This section gives an overview of all files that are being used during processing, how to implement and when to use them.

### fastq

This is the standard output format for sequencing data.
It contains all (raw) sequencing reads in random order including a quality string per basepair.
Each read has four lines:

1. Header: This contains some basic information from the sequencing machine and a unique identifier number.
2. Sequence: This is the actual nucleotide sequence.
3. Dummy: This is typically a '+' and is there to separate the sequence line from the quality line.
4. Quality score: This indicates the quality of each basepair in the sequence line (each symbol in this line belongs the nucleotide at the same position in the sequence line). The sequence and this quality line should always have the same length. The quality score is given in terms of phred scores (see below).

The quality line is typically given as a phred score.
There are two versions, base33 and base64, but the base64 is outdated and hardly used anymore.
In both versions the quality score is determined by Q = -10*log10(P) where P is the error probability (0 < P < 1).
A Q-score of 0 (i.e. an error probability of P=1) is defined by ascii symbol 33 ('!') for base33 and by ascii symbol 64 ('@') for base64.
A Q-score of 1 (p=0.79) is then given by ascii 34 (' " ') (for base33) etcetera.
For a full table of ascii symbols and probability scores, see the appedices of this document [PHRED tables (base33)](#phred-tables-base33) and [PHRED tables (base64)](#phred-tables-base64).

The nucleotide sequence typically only contains the four nucleotide letters (A, T, C and G), but when a nucleotide was not accurately determined (i.e. having a error probability higher than a certain threshold), the nucleotide is sometimes converted to the letter N, indicating that this nucleotide was not successfully sequenced.

Fastq files tend to be large in size (depending on how many reads are sequenced, but >10Gb is normal).
Therefore these files are typically compressed in gzip format (.fastq.gz).
The pipeline can handle gzipped files by itself, so there is no need to convert it manually.

Example:

> `@NB501605:544:HLHLMBGXF:1:11101:9938:1050 1:N:0:TGCAGCTA`  
> `TGTCAACGGTTTAGTGTTTTCTTACCCAATTGTAGAGACTATCCACAAGGACAATATTTGTGACTTATGTTATGCG`  
> `+`  
> `AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE`  
> `@NB501605:544:HLHLMBGXF:1:11101:2258:1051 1:N:0:TACAGCTA`  
> `TGAGGCACCTATCTCAGCGATCGTATCGGTTTTCGATTACCGTATTTATCCCGTTCGTTTTCGTTGCCGCTATTT`  
> `+`  
> `AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEEEEEE<EEEAEEEEEEE/E/E/EA///`  
> `@NB501605:544:HLHLMBGXF:1:11101:26723:1052 1:N:0:TGCAGCTA`  
> `TGTCAACGGTTTAGTGTTTTCTTACCCAATTGTAGAGACTATCCACAAGGACAATATTTGTGACTTATGTTATGCG`  
> `+`  
> `AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE`

### sam & bam

### wig

### bed

### pergene.txt & peressential.txt

### pergene_insertions.txt & peressential_insertions.txt

## Software - Processing

### satay.sh

- **Main tasks**

Briefly explain what this thing does.

- **Dependencies**

All scripts and files where this thing is depending on.

- **How does it work**

A more detailed explanation how the thing works and optionally different functions within the thing.

- **How to use**

What to input, how to run and arguments and settings, what to expect from output, how to change things, what to do in case of error.

- **Output files**

More detailed explanation about output if required.

- **Notes**

Some extra note to be aware of.

## Software - analysis

### python scripts

#### strip_redundant_insertions.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### genomicfeatures_dataframe.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### scatterplot_genes.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### transposonread_profileplot.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### transposonread_profileplot_genome.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### volcanoplot.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### create_essentialgenes_list.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

#### split_wigfiles.py

- **Main tasks**

- **Dependencies**

- **How does it work**

- **How to use**

- **Output files**

- **Notes**

### python modules

#### chromosome_and_gene_positions.py

#### chromosome_names_in_files.py

#### dataframe_from_pergene.py

#### essential_genes_names.py

#### gene_length.py

#### gene_names.py

#### gene_tn_insertions.py

#### insertions_count.py

#### mapped_reads.py

#### read_sgdfeatures.py

#### statistics_perchromosome.py

### Other tools

#### IGV

#### genome browser

## Outlook

## How to use the Linux desktop

general layout

Log in codes

sftp for drives

command line

## Appendices

### PHRED tables (base33)

| ASCII symbol | ASCII value | Q  | P\_err   |
| ------------ | ----------- | -- | -------- |
| !            | 33          | 0  | 1,000000 |
| "            | 34          | 1  | 0,794328 |
| #            | 35          | 2  | 0,630957 |
| $            | 36          | 3  | 0,501187 |
| %            | 37          | 4  | 0,398107 |
| &            | 38          | 5  | 0,316228 |
|  '           | 39          | 6  | 0,251189 |
| (            | 40          | 7  | 0,199526 |
| )            | 41          | 8  | 0,158489 |
| \*           | 42          | 9  | 0,125893 |
| +            | 43          | 10 | 0,100000 |
| ,            | 44          | 11 | 0,079433 |
| \-           | 45          | 12 | 0,063096 |
| .            | 46          | 13 | 0,050119 |
| /            | 47          | 14 | 0,039811 |
| 0            | 48          | 15 | 0,031623 |
| 1            | 49          | 16 | 0,025119 |
| 2            | 50          | 17 | 0,019953 |
| 3            | 51          | 18 | 0,015849 |
| 4            | 52          | 19 | 0,012589 |
| 5            | 53          | 20 | 0,010000 |
| 6            | 54          | 21 | 0,007943 |
| 7            | 55          | 22 | 0,006310 |
| 8            | 56          | 23 | 0,005012 |
| 9            | 57          | 24 | 0,003981 |
| :            | 58          | 25 | 0,003162 |
| ;            | 59          | 26 | 0,002512 |
| <            | 60          | 27 | 0,001995 |
| \=           | 61          | 28 | 0,001585 |
| \>           | 62          | 29 | 0,001259 |
| ?            | 63          | 30 | 0,001000 |
| @            | 64          | 31 | 0,000794 |
| A            | 65          | 32 | 0,000631 |
| B            | 66          | 33 | 0,000501 |
| C            | 67          | 34 | 0,000398 |
| D            | 68          | 35 | 0,000316 |
| E            | 69          | 36 | 0,000251 |
| F            | 70          | 37 | 0,000200 |
| G            | 71          | 38 | 0,000158 |
| H            | 72          | 39 | 0,000126 |
| I            | 73          | 40 | 0,000100 |
| J            | 74          | 41 | 0,000079 |
| K            | 75          | 42 | 0,000063 |
| L            | 76          | 43 | 0,000050 |
| M            | 77          | 44 | 0,000040 |
| N            | 78          | 45 | 0,000032 |
| O            | 79          | 46 | 0,000025 |
| P            | 80          | 47 | 0,000020 |
| Q            | 81          | 48 | 0,000016 |
| R            | 82          | 49 | 0,000013 |
| S            | 83          | 50 | 0,000010 |

### PHRED table (base64)

| ASCII symbol | ASCII value | Q  | P\_err   |
| ------------ | ----------- | -- | -------- |
| @            | 64          | 0  | 1,000000 |
| A            | 65          | 1  | 0,794328 |
| B            | 66          | 2  | 0,630957 |
| C            | 67          | 3  | 0,501187 |
| D            | 68          | 4  | 0,398107 |
| E            | 69          | 5  | 0,316228 |
| F            | 70          | 6  | 0,251189 |
| G            | 71          | 7  | 0,199526 |
| H            | 72          | 8  | 0,158489 |
| I            | 73          | 9  | 0,125893 |
| J            | 74          | 10 | 0,100000 |
| K            | 75          | 11 | 0,079433 |
| L            | 76          | 12 | 0,063096 |
| M            | 77          | 13 | 0,050119 |
| N            | 78          | 14 | 0,039811 |
| O            | 79          | 15 | 0,031623 |
| P            | 80          | 16 | 0,025119 |
| Q            | 81          | 17 | 0,019953 |
| R            | 82          | 18 | 0,015849 |
| S            | 83          | 19 | 0,012589 |
| T            | 84          | 20 | 0,010000 |
| U            | 85          | 21 | 0,007943 |
| V            | 86          | 22 | 0,006310 |
| W            | 87          | 23 | 0,005012 |
| X            | 88          | 24 | 0,003981 |
| Y            | 89          | 25 | 0,003162 |
| Z            | 90          | 26 | 0,002512 |
| [            | 91          | 27 | 0,001995 |
| \            | 92          | 28 | 0,001585 |
| ]            | 93          | 29 | 0,001259 |
| ^            | 94          | 30 | 0,001000 |
| _            | 95          | 31 | 0,000794 |
| `            | 96          | 32 | 0,000631 |
| a            | 97          | 33 | 0,000501 |
| b            | 98          | 34 | 0,000398 |
| c            | 99          | 35 | 0,000316 |
| d            | 100         | 36 | 0,000251 |
| e            | 101         | 37 | 0,000200 |
| f            | 102         | 38 | 0,000158 |
| g            | 103         | 39 | 0,000126 |
| h            | 104         | 40 | 0,000100 |
| i            | 105         | 41 | 0,000079 |
| j            | 106         | 42 | 0,000063 |
| k            | 107         | 43 | 0,000050 |
| l            | 108         | 44 | 0,000040 |
| m            | 109         | 45 | 0,000032 |
| n            | 110         | 46 | 0,000025 |
| o            | 111         | 47 | 0,000020 |
| p            | 112         | 48 | 0,000016 |
| q            | 113         | 49 | 0,000013 |
| r            | 114         | 50 | 0,000010 |

> *"I want to be a healer, and love all things that grow and are not barren"* - J.R.R. Tolkien
