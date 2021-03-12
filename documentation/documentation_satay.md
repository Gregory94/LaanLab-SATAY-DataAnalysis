# Documentation for processing SATAY <!-- omit in toc -->

- [Introduction](#introduction)
- [File types](#file-types)
  - [fastq](#fastq)
  - [sam & bam](#sam--bam)
  - [bed](#bed)
  - [wig](#wig)
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
    - [samflag.py](#samflagpy)
  - [Other tools](#other-tools)
    - [IGV](#igv)
    - [genome browser](#genome-browser)
- [Outlook](#outlook)
- [How to use the Linux desktop](#how-to-use-the-linux-desktop)
- [Appendices](#appendices)
  - [PHRED table (base33)](#phred-table-base33)
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
This section gives an overview of all those files, how to implement them and when to use which file.

### fastq

This is the standard output format for sequencing data.
It contains all (raw) sequencing reads in random order including a quality string per basepair.
Each read has four lines:

1. Header: Contains some basic information from the sequencing machine and a unique identifier number.
2. Sequence: The actual nucleotide sequence.
3. Dummy: Typically a '+' and is there to separate the sequence line from the quality line.
4. Quality score: Indicates the quality of each basepair in the sequence line (each symbol in this line belongs to the nucleotide at the same position in the sequence line). The sequence and this quality line should always have the same length.

The quality line is given as a phred score.
There are two versions, base33 and base64, but the base64 is outdated and hardly used anymore.
In both versions the quality score is determined by Q = -10*log10(P) where P is the error probability determined during sequencing (0 < P < 1).
A Q-score of 0 (i.e. an error probability of P=1) is defined by ascii symbol 33 ('!') for base33 and by ascii symbol 64 ('@') for base64.
A Q-score of 1 (p=0.79) is then given by ascii 34 (' " ') (for base33) etcetera.
For a full table of ascii symbols and probability scores, see the appendices of this document [PHRED table (base33)](#phred-table-base33) and [PHRED table (base64)](#phred-table-base64).

The nucleotide sequence typically only contains the four nucleotide letters (A, T, C and G), but when a nucleotide was not accurately determined (i.e. having a error probability higher than a certain threshold), the nucleotide is sometimes converted to the letter N, indicating that this nucleotide was not successfully sequenced.

Fastq files tend to be large in size (depending on how many reads are sequenced, but >10Gb is normal).
Therefore these files are typically compressed in gzip format (.fastq.gz).
The pipeline can handle gzipped files by itself, so there is no need to convert it manually.

Example fastq file:

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

When the reads are aligned to a reference genome, the resulting file is a Sequence Alignment Mapping (sam) file.
Every read is one line in the file and consists of at least 11 tab delimited fields that are always in the same order:

1. Name of the read. This is unique for each read, but can occur multiple times in the file when a read is split up over multiple alignments at different locations.
2. Flag. Indicating properties of the read about how it was mapped (see below for more information).
3. Chromosome name to which the read was mapped.
4. Leftmost position in the chromosome to which the read was mapped.
5. Mapping quality in terms of Q-score as explained in the [fastq](#fastq) section.
6. CIGAR string. Describing which nucleotides were mapped, where insertions and deletions are and where mismatches occurs. For example, `43M1I10M3D18M` means that the first 43 nucleotides match with the reference genome, the next 1 nucleotide exists in the read but not in the reference genome (insertion), then 10 matches, then 3 nucleotides that do not exist in the read but do exist in the reference genome (deletions) and finally 18 matches. For more information see [this website](https://www.drive5.com/usearch/manual/cigar.html).
7. Reference name of the mate read (when using paired end datafiles). If no mate was mapped (e.g. in case of single end data or if it was not possible to map the mate read) this is typically set to `*`.
8. Position of the mate read. If no mate was mapped this is typically set to `0`.
9. Template length. Length of a group (i.e. mate reads or reads that are split up over multiple alignments) from the left most base position to the right most base position.
10. Nucleotide sequence.
11. Phred score of the sequence (see [fastq](#fastq) section).

Depending on the software, the sam file typically starts with a few header lines containing information regarding the alignment.
For example for BWA MEM (which is used in the pipeline), the sam file start with `@SQ` lines that shows information about the names and lengths for the different chromosome and `@PG` shows the user options that were set regarding the alignment software.
Note that these lines might be different when using a different alignment software.
Also, there is a whole list of optional fields that can be added to each read after the first 11 required fields.
For more information, see [wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)).

The flag in the sam files is a way of representing a list of properties for a read as a single integer.
There is a defined list of properties in a fixed order:

1. read paired
2. read mapped in proper pair
3. read unmapped
4. mate unmapped
5. read reverse strand
6. mate reverse strand
7. first in pair
8. second in pair
9. not primary alignment
10. read fails platform/vendor quality checks
11. read is PCR or optical duplicate
12. supplementary alignment

To determine the flag integer, a 12-bit binary number is created with zeros for the properties that are not true for a read and ones for those properties that are true for that read.
This 12-bit binary number is then converted to a decimal integer.
Note that the binary number should be read from right to left.
For example, FLAG=81 corresponds to the 12-bit binary 000001010001 which indicates the properties: 'read paired', 'read reverse strand' and 'first in pair'.
Decoding of sam flags can be done using [this website](http://broadinstitute.github.io/picard/explain-flags.html) or using [samflag.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/python_modules/samflag.py "LaanLab-SATAY_DataAnalysis.samflag.py").

Example sam file (note that the last read was not mapped):

>`NB501605:544:HLHLMBGXF:1:11101:25386:1198     2064    ref|NC_001136|  362539  0       42H34M  *       0       0       GATCACTTCTTACGCTGGGTATATGAGTCGTAAT      EEEAEEEEEEAEEEAEEAEEEEEEEEEEEAAAAA      NM:i:0  MD:Z:34 AS:i:34 XS:i:0  SA:Z:ref|NC_001144|,461555,+,30S46M,0,0;`
>
>`NB501605:544:HLHLMBGXF:1:11101:20462:1198       16      ref|NC_001137|  576415  0       75M     *       0       0       CTGTACATGCTGATGGTAGCGGTTCACAAAGAGCTGGATAGTGATGATGTTCCAGACGGTAGATTTGATATATTA     EEEAEEEEEEEEEEEEEAEEAE/EEEEEEEEEEEAEEEEEE/EEEEEAEAEEEEEEEEEEEEEEEEEEEEAAAAA     NM:i:1  MD:Z:41C33      AS:i:41 XS:i:41`
>
>`NB501605:544:HLHLMBGXF:1:11101:15826:1199       4       *       0       0       *       *       0       0       ACAATATTTGTGACTTATGTTATGCG      EEEEEEEEEEE6EEEEEEEEEEEEEE      AS:i:0  XS:i:0`

Sam files tend to be large in size (tens of Gb is normal).
Therefore the sam files are typically stored as compressed binary files called bam files.
Almost all downstream analysis tools that need the alignment information accept bam files as input.
Therefore the sam files are mostly deleted after the bam file is created.
When a sam file is needed, it can always be recreated from the bam file, for example using `samtools` using the command `samtools view -h -o out.sam in.bam`.

### bed

A bed file is one of the outputs from the transposon mapping pipeline.
It is a standard format for storing read insertion locations and the corresponding read counts.
The file consists of a single header, typically something similar to `track name=[file_name] userscore=1`.
Every row corresponds to one insertion and has (in case of the satay analysis) the following space delimited columns:

1. chromosome (e.g. `chrI` or `chrref|NC_001133|`)
2. start position of insertion
3. end position of insertion (in case of satay-analysis, this is always start position + 1)
4. dummy column (this information is not present for satay analysis, but must be there to satisfy the bed format)
5. number of reads at that insertion location

In case of processing with `transposonmapping.py` (final step in processing pipeline) or [the matlab code from the kornmann-lab](https://sites.google.com/site/satayusers/complete-protocol/bioinformatics-analysis/matlab-script), the number of reads are given according to `(reads*20)+100`, for example 2 reads are stored as 140.

The bed file can be used for many downstream analysis tools, for example [genome_browser](http://genome-euro.ucsc.edu/index.html).

Sometimes it might occur that insertions are stored outside the chromosome (i.e. the insertion position is bigger than the length of that chromosome).
Also, reference genomes sometimes do not have the different chromosomes stored as roman numerals (for example `chrI`, `chrII`, etc.) but rather use different names (this originates from the chromosome names used in the reference genome).
These things can confuse some analysis tools, for example the [genome_browser](http://genome-euro.ucsc.edu/index.html).
To solve this, the python function [strip_redundant_insertions.py](#strip_redundant_insertionspy) is created.
This creates a _clean.bed file where the insertions outside the chromosome are removed and all the chromosome names are stored with their roman numerals.
See [strip_redundant_insertions.py](#strip_redundant_insertionspy) for more information.

Example bed file:

> `track name=leila_wt_techrep_ab useScore=1`  
> `chrI 86 87 . 140`  
> `chrI 89 90 . 140`  
> `chrI 100 101 . 3820`  
> `chrI 111 112 . 9480`

### wig

A wiggle (wig) file is another output from the transposon mapping pipeline.
It stores similar information as the bed file, but in a different format.
This file has a header typically in the form of `track type=wiggle_0 maxheightPixels=60 name=[file_name]`.
Each chromosome starts with the line `variablestep chrom=chr[chromosome]` where `[chromosome]` is replaced by a chromosome name, e.g. `I` or `ref|NC_001133|`.
After a variablestep line, every row corresponds with an insertion in two space delimited columns:

1. insertion position
2. number of reads

In the wig file, the read count represent the actual count (unlike the bed file).

There is one difference between the bed and the wig file.
In the bed file the insertions at the same position but with a different orientation are stored as individual insertions.
In the wig file these insertions are represented as a single insertion and the corresponding read counts are added up.

Similar to the bed file, also in the wig insertions might occur that have an insertion position that is bigger then the length of the chromosome.
This can be solved with the [same python script](#strip_redundant_insertionspy) as the bed file.
The insertions that have a position outside the chromosome are removed and the chromosome names are represented as a roman numeral.

Example wig file:

> `track type=wiggle_0 maxheightPixels=60 name=WT_merged-techrep-a_techrep-b_trimmed.sorted.bam`  
> `variablestep chrom=chrI`  
> `86 2`  
> `89 2`  
> `100 186`  
> `111 469`

### pergene.txt & peressential.txt

A pergene.txt and peressential.txt file are yet another outputs from the transposon mapping pipeline.
Where bed and wig files store *all* insertions throughout the genome, these files only store the insertions in each gene or each essential gene, respectively.
Essential genes are the annotated essential genes as stated by SGD for wild type cells.
The genes are taken from the [Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Yeast_Protein_Names.txt) file, which is downloaded from [uniprot](https://www.uniprot.org/docs/yeast).
The positions of each gene are determined by a [gff3 file](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Saccharomyces_cerevisiae.R64-1-1.99.gff3) downloaded from SGD.
Essential genes are defined in [Cerevisiae_AllEssentialGenes.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_AllEssentialGenes_List.txt).

The pergene.txt and the peressential.txt have the same format.
This consists of a header and then each row contains three tab delimited columns:

1. gene name
2. total number of insertions within the gene
3. sum of all reads of those insertions

The reads are the actual read counts.

Note that when comparing files that include gene names there might be differences in the gene naming.
Genes have multiple names, e.g. systematic names like 'YBR200W' or standard names like 'BEM1' which can have aliases such as 'SRO1'.
The above three names all refer to the same gene.
The [Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Yeast_Protein_Names.txt) file can be used to search for aliases when comparing gene names in different files, or the [genomicfeatures_dataframe.py](#genomicfeatures_dataframepy) python script can be used which creates a pandas dataframe that includes the different gene names (this python script itself makes also use of the Yeast_Protein_Names.txt file).

Example of pergene.txt file:

> `Gene name Number of transposons per gene Number of reads per gene`  
> `YAL069W 34 1819`  
> `YAL068W-A 10 599`  
> `PAU8 26 1133`  
> `YAL067W-A 12 319`

### pergene_insertions.txt & peressential_insertions.txt

The final two files that are created by the transposon mapping pipeline are the pergene_insertions.txt and the peressential_insertions.txt.
The files have a similar format as the pergene.txt file, but are more extensive in terms of the information per gene.
The information is taken from [Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Yeast_Protein_Names.txt), the [gff3 file](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Saccharomyces_cerevisiae.R64-1-1.99.gff3) and [Cerevisiae_AllEssentialGenes.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/Data_Files/Cerevisiae_AllEssentialGenes_List.txt), similar as the pergene.txt files.

Both the pergene_insertions.txt and the peressential_insertions.txt files have a header and then each row contains six tab delimited columns:

1. gene name
2. chromosome where the gene is located
3. start position of the gene
4. end position of the gene
5. list of all insertions within the gene
6. list of all read counts in the same order as the insertion list

This file can be useful when not only the number of insertions is important, but also the distribution of the insertions within the genes.

Example of pergene_insertions.txt file:

> `Essential gene name chromosome Start location End location Insertion locations Reads per insertion location`  
> `EFB1 I 142174 143160 [142325, 142886] [1, 1]`  
> `MAK16 I 100225 101145 100229, 100407, 100422, 100791, 101022, 101129] [4, 1, 5, 1, 1, 1]`  
> `PRE7 II 141247 141972 [141262, 141736, 141742, 141895] [1, 1, 1, 1]`  
> `RPL32 II 45978 46370 [46011, 46142, 46240] [1, 3, 1]`

## Software - Processing

This section discusses the main pipeline for processing satay datasets, from the raw fastq files of the sequencing output to the bed, wig and pergene text files which can be directly used for analysis.
See the image below for a schematic overview of the pipeline with the used software tools between brackets and the file type after each processing step on the left.
Everything that is shown in green can be automatically performed in a single workflow as will be discussed in [satay section](#sataysh) below.

<img src="C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\documentation\media\satayprocessingpipeline.png" alt="satay_processing_pipeline" width="500"/>

### satay.sh

*Note that the workflow only runs in Linux (or Mac, but this is untested) as some of its dependencies are specifically designed for Unix machines.*

The main program for the data processing is called [satay.sh](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/satay.sh) and is a bash script which automatically calls all required software tools.

- **Dependencies**

Below is a list of dependencies.
Note the folder structure for the python scripts and modules and the data files.

1. Zenity (by default installed in Linux Ubuntu 20.10)
2. [YAD](http://manpages.ubuntu.com/manpages/groovy/man1/yad.1.html)
3. [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/)
5. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
6. [BWA](http://bio-bwa.sourceforge.net/bwa.shtml)
7. [SAMTools](http://www.htslib.org/)
8. [Sambamba](https://lomereiter.github.io/sambamba/)
9. Python > v3.7
    a. numpy
    b. [pysam](https://pysam.readthedocs.io/en/latest/index.html)
    c. timeit
    d. python_scripts/[transposonmapping_satay.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/transposonmapping_satay.py)
    e. python_scripts/python_modules/[chromosome_and_gene_positions.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/chromosome_and_gene_positions.py)
    f. python_scripts/python_modules/[gene_names.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/gene_names.py)
    g. python_scripts/python_modules/[samflag.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/samflag.py)
10. resources/Saccharomyces_cerevisiae.R64-1-1.99.gff3
11. resources/Cerevisiae_AllEssentialGenes_List.txt
12. resources/Yeast_Protein_Names.txt

- **Main tasks**

Briefly explain what this thing does.

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
  
Difference in VariableStep in variablestep

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

#### samflag.py

### Other tools

#### IGV

#### genome browser

## Outlook

## How to use the Linux desktop

general layout

Log in codes

sftp for drives and disabling after 10 minutes

command line

updating

## Appendices

### PHRED table (base33)

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
