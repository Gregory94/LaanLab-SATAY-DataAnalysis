# Documentation for processing SATAY <!-- omit in toc -->

- [Introduction](#introduction)
- [File types](#file-types)
  - [fastq](#fastq)
  - [sam, bam](#sam-bam)
  - [bed](#bed)
  - [wig](#wig)
  - [pergene.txt, peressential.txt](#pergenetxt-peressentialtxt)
  - [pergene_insertions.txt, peressential_insertions.txt](#pergene_insertionstxt-peressential_insertionstxt)
- [Software - Processing](#software---processing)
  - [satay.sh](#sataysh)
    - [Dependencies](#dependencies)
    - [Input, Output](#input-output)
    - [How to use](#how-to-use)
    - [Tutorial](#tutorial)
    - [How does it work](#how-does-it-work)
    - [Notes](#notes)
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
Basically all fastq files that are being created on modern sequencing machines use the base33 system.

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

### sam, bam

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
Almost all downstream analysis tools (at least all tools discussed in this document) that need the alignment information accept bam files as input.
Therefore the sam files are mostly deleted after the bam file is created.
When a sam file is needed, it can always be recreated from the bam file, for example using `SAMTools` using the command `samtools view -h -o out.sam in.bam`.
The bam file can be sorted (creating a .sorted.bam file) where the reads are typically ordered depending on their position in the genome.
This usually also comes with an index file (a .sorted.bam.bai file) which stores some information where for example the different chromosomes start within the bam file and where specific often occuring sequences are.
Many downstream analysis tools require this file to be able to efficiently search through the bam file.

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

In the wig file, the read count represent the actual count (unlike the bed file where an equation is used to transform the numbers).

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

### pergene.txt, peressential.txt

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
To suppress noise, the insertion with the highest read count in a gene is removed from that gene.
Therefore, it might occur that a gene has 1 insertion, but 0 reads.

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

### pergene_insertions.txt, peressential_insertions.txt

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
Similarily as the [pergene.txt and peresential.txt file](#pergenetxt-peressentialtxt), to suppress noise the insertion with the highest read count in a gene is removed from that gene.

Example of pergene_insertions.txt file:

> `Essential gene name chromosome Start location End location Insertion locations Reads per insertion location`  
> `EFB1 I 142174 143160 [142325, 142886] [1, 1]`  
> `MAK16 I 100225 101145 100229, 100407, 100422, 100791, 101022, 101129] [4, 1, 5, 1, 1, 1]`  
> `PRE7 II 141247 141972 [141262, 141736, 141742, 141895] [1, 1, 1, 1]`  
> `RPL32 II 45978 46370 [46011, 46142, 46240] [1, 3, 1]`

## Software - Processing

This section discusses the main pipeline for processing satay datasets, from the raw fastq files of the sequencing output to the bed, wig and pergene text files which can be directly used for analysis.
See the image below for a schematic overview of the pipeline with the used software tools between brackets and the file type after each processing step on the left.
Everything that is shown in green can be automatically performed in a single workflow as will be discussed in the [satay section](#sataysh) below.

<img src="C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\documentation\media\satayprocessingpipeline.png" alt="satay_processing_pipeline" width="500"/>

### satay.sh

*Note that the workflow only runs in Linux (or Mac, but this is untested) as some of its dependencies are specifically designed for Unix machines.*

The main program for the data processing is called [satay.sh](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/satay.sh) and is a bash script which automatically calls all required software tools.

#### Dependencies

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
   1. numpy
   2. [pysam](https://pysam.readthedocs.io/en/latest/index.html)
   3. timeit
   4. python_scripts/[transposonmapping_satay.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/transposonmapping_satay.py)
   5. python_scripts/python_modules/[chromosome_and_gene_positions.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/chromosome_and_gene_positions.py)
   6. python_scripts/python_modules/[gene_names.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/gene_names.py)
   7. python_scripts/python_modules/[samflag.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/python_modules/samflag.py)
10. data_files/[Saccharomyces_cerevisiae.R64-1-1.99.gff3](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/data_files/Saccharomyces_cerevisiae.R64-1-1.99.gff3)
11. data_files/[Cerevisiae_AllEssentialGenes_List.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/data_files/Cerevisiae_AllEssentialGenes_List.txt)
12. data_files/[Yeast_Protein_Names.txt](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/data_files/Yeast_Protein_Names.txt)

#### Input, Output

The bash script satay.sh accepts raw fastq files, either compressed (gzip) or uncompressed, with one of the following extensions: .fastq, .fq, fastq.gz, fq.gz.
The fastq files can be either paired-end or single-end.

It is assumed that each data file contains one sample, so in case multiple samples are sequenced together, it might be necessary to demultiplex the fastq files.
This is not integrated in the pipeline and should be performed before the pipeline is started.

When the processing was completed successfully, a number of output folders are created at the same location where the input fastq file is stored. The creation of some files and folders are depending on which options are selected in satay.sh (see below for more information about the options).
The output may include the following folders and files:

1. `align_out/`
   1. **.sam** (depending on whether the option for deleting sam files was selected)
   2. **.bam**
   3. **.sorted.bam** (depending on whether the option for `Sort and index bam files` is selected); The reads are sorted which speeds up downstream processing.
   4. **.sorted.bam.bai** (depending on whether the option for `Sort and index bam files` is selected);
   5. **flagstatreport.txt** (depending on whether `Create flagstat report is selected` is selected); Stores some basic information about the alignment.
   6. **.bed** (depending on whether `Transposon mapping` is selected)
   7. **.wig** (depending on whether `Transposon mapping` is selected)
   8. **pergene.txt** (depending on whether `Transposon mapping` is selected)
   9. **peressential.txt** (depending on whether `Transposon mapping` is selected)
   10. **pergene_insertions.txt** (depending on whether `Transposon mapping` is selected)
   11. **peressential_insertions.txt** (depending on whether `Transposon mapping` is selected)
2. `trimm_out/` (depending on whether the option for trimming is selected)
   1. **trimmed.fastq** (depending on whether trimming is selected)
3. `fastqc/` (depending on whether the option for quality checking is selected)
   1. **.html** (depending on whether Quality checking is selected); contains an a number of figures showing the quality of the data. This can be created for both raw data and/or trimmed data.
   2. **zipped folder** (depending on whether Quality checking is selected); contains the data for creating the figures in the .html file.

#### How to use

The workflow satay.sh only runs in Linux and is designed to be used as a commandline tool (see [How to use the Linux desktop](#how-to-use-the-linux-desktop)).
For a step-by-step guide, see the [Tutorial](#tutorial).

In the commandline, navigate to the software folder (`cd /home/laanlab/Documents/satay/software/`).
The workflow can be started either using commandline arguments or by using a graphical user interface (GUI).
Using the GUI is the most userfriendly approach.
Access the help text using `bash satay.sh --help` or `bash satay.sh -h` and check the current version with `bash satay.sh -v`.
The help text explains the different commandline arguments that can be set and briefly how to use the workflow.
The commandline arguments are the same options that are set when using the GUI, so here only the GUI is explained.

Start the GUI with `bash satay.sh` without any arguments.
This opens a window where the data file(s) can be selected.
Navigate to the folder containing the data file(s) and select one or two files (select only two files when using paired-end data that is not interleaved by holding the ctrl button and clicking the files).
Use the drop down window at the bottom right to select the right extension.

<img src="C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\documentation\media\satay_fileselectionwindow.png" alt="satay.sh_file_selection_window" width=700>

After pressing `OK` a second window opens that shows all the settings that can be changed with some default values.
The first two lines indicate the file(s) that were selected in the previous window (if only one file was selected, the secondary reads file is set to `none`).
The third line is a drop-down menu to set whether to process the data as single-end (default) or paired-end reads.
The fourth line sets which trimming tool to use.
There are two trimming tools installed, BBDuk and Trimmomatic, which can be selected in this drop-down menu (default is BBDuks).
If the trimming should be skipped (and therefore the alignment is performed with the raw input data), select `donottrim` in this menu.
The fifth line sets the trimming arguments.
When trimming is skipped, what is put in this line is ignored.
Otherwise, use the arguments allowed by the selected trimming software.

**Trimming settings BBDuk**: For BBDuk, a typical argument line looks like this: `ktrim=l k=15 mink=10 hdist=1 qtrim=r trimq=10 minlen=30`.
The first four options are for trimming unwanted (adapter) sequences which in BBDuk is done using k-mers (e.g. when k=3, all 3-mers for ATTGCAAT are ATT, TTG, TGC, GCA, CAA, AAT).
Any unwanted sequences that need to be trimmed (e.g. adapter, primer and transposon sequences) should be entered in the adapters file.
Clicking the `Open adapters file` opens a text file in which these unwanted sequences can be placed in fasta format.
Fasta format is similar to the fastq format, but without the last two lines (i.e. the dummy line and the quality line).
The header line of each read should start with a `>` and can contain any text (preferably no special symbols to be sure) and the next line should contain the unwanted sequence, for example:

> `> header line first unwanted sequence`  
> `GATC`  
> `> header line second unwanted seqence`  
> `CATG`

When all unwanted sequences are set, save the text file and close it.
The value for k is set by the `k` parameter and should not be higher than the length of the shortest unwanted sequence (e.g. if the shortest sequence that needs to be trimmed has length 10, then k < 10).
The higher this number, the smaller the chances are that false positives are found when trimming.
This method, however, might skip over the last few basepairs of a read when the length of a read is not a multiple of k.
To solve for this, the `mink` option can be set which allows for shorter k-mers at the end of a read.
Therefore, `mink` < `k`.
The `ktrim` option set to which side to trim when encountering an unwanted sequence (right (`r`) or left (`l`)).
The number set with `hdist` determines the number of mistakes allowed when matching the unwanted sequences with the reads.
Next are some options for quality trimming of the reads.
The `trimq` option sets the minimum Q-score that basepairs are allowed to have on average.
When the average quality of either first few or last few basepairs is lower than this threshold, those basepairs will be removed.
Whether the basepairs on the left or on the right have to be trimmed is determined by the `qtrim` parameter (either `r` or `l` or `rl` to trim everything in that read).
Finally the `minlen` parameter checks for the minimum length of a read after trimming.
When this is smaller then this threshold, the read is completely removed.
Note that the trimming can have a serious influence on the alignment and choosing good options can be important.
For more information about these and many other settings check the [BBDuk manual](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

**Trimming settings Trimmomatic**: Alternatively, instead of using BBDuk for the trimming, Trimmomatic can be used that requires different input arguments.
A typical input for Trimmomatic can look like this (note the capital letters): `ILLUMINACLIP:1:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:5:15 MINLEN:10`.
The `ILLUMINACLIP` argument performs the trimming of the unwanted sequences from the adapters file.
Use the same file as for BBDuk which can be accessed using the `Open adapters file` button at the bottom of the window and also enter the reads in the same fasta format.
The first number of this argument indicates how many mismatches are allowed.
The second number is mainly for paired-end reads and determines how accurate a match should be when considering a read pair (even though this is only for paired-end reads, it is required to enter this value, but it is not important for single-end reads).
The third value sets the accuracy of a match between an unwanted sequence and the read.
The `SLIDINGWINDOW` is a quality control and uses a sliding window approach to determine the average quality of the reads within the window (length in basepairs is given by the first value) and trim the read when the quality drops below a threshold given by the second value.
The `LEADING` and `TRAILING` arguments check the quality of the basepairs at the beginning or end of the read, respectively.
If the quality is below the threshold given by the value, that basepair is trimmed.
This is repeated for all basepairs from the start of the read or at the end of the read until a basepair is found that has a quality above the threshold.
The `MINLEN` argument sets the minimum length of the read after trimming.
If the length is below the threshold, the entire read is removed.
For more information about these and many other settings check the [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

**Alignment settings**: The alignment settings line takes the arguments for the alignment software.
The default value is `-v 2` which reduces the verbose level and prevents a whole lot of text to be printed during alignment, but doesn't do anything with the way the aligner works.
The aligner determines for each read the position in the genome where it has the highest mapping score.
The read is aligned to a location and then this score is determined by giving points for each letter in the read that match with the reference genome, but adding penalty points for evey time a letter in the read doesn't match with the genome (e.g. because of mismatches, extra letters in the read or missing letters in the read relative to the reference genome).
This score is determined for each location where the read is aligned and the read is mapped to the location with the highest score.
When the read cannot be unambiguously aligned (e.g. because all the mapping scores are too low or when there are mulitple locations that have the same highest score), the reads are aligned to a random location.
The most important parameters that can be set alters the way the mapping score is influenced by matches, mismatches, insertions, deletions etc.
Useful parameters to consider are: `-A`, `-B`, `-O`, `-E` and `-L`.
When using paired end data it is also good to think about using the parameter `-U` to change the penalty for unpaired read pairs.
For explanation about the parameters and their default values, see the [BWA MEM manual](http://bio-bwa.sourceforge.net/bwa.shtml).

Next are check boxes that can be set to turn on or off certain parts of the workflow.
The first two checkboxes turn on or off quality checking before and after trimming.
When `Quality check interrupt` is set, the program pauses after the raw quality checking (before the trimming) and asking if the user want to continue the processing (press `n` to stop and `y` to continue).
Stopping the program allows the user to check the quality report of the raw data.
After checking the quality report, the program can be restarted by typing `bash satay.sh` after which the program remembers the settings that were previously set (it will skip the file selection window).
Changes can be made for the trimming and alignment tools and press `OK` to continue.
This time the raw quality report will be skipped.

The next options determine if the sam file needs to be deleted after processing, whether the bam file needs to be sorted and indexed, if the [transposon mapping](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/transposonmapping_satay.py) needs to be performed and if a flagstat report needs to be created (the quality report of the alignment).

<img src="C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\documentation\media\satay_settingswindow.png" alt="satay.sh_settings_window" width=700>

#### Tutorial

**Note before use**: When down-/uploading files to the N-drive or M-drive on the Linux Desktop, the drives will disconnect automatically after 10 minutes of inactivity.
Since the files are typically large the down-/uploading can therefore take more then 10 minutes to complete.
To prevent disconnection while down-/uploading you should click back and forth to some folders once every few minutes in the drives to reset the timer.
See also [How to use the Linux desktop](#how-to-use-the-linux-desktop).

- [ ] Log in to the Linux desktop (log in credentials can be found at `N:\tnw\BN\LL\Shared\LinuxMachines\LinuxDesktop_LoginCredentials.txt`).
  - [ ] When using remote access, install the free version of [Teamviewer](https://www.teamviewer.com/nl/) on your computer and connect with the computer using the credentials found in the file mentioned above.
- [ ] If your data file is on the N-drive, it should be copied to the Linux desktop for the processing to work smoothly (if you have no data set, there is a test fastq file already on the Linux desktop, so you can skip the following steps):
  - [ ] For this, go to `Files` (located in the left bar) and click `Other Locations` on the left menu in Files.
  - [ ] On the bottom at `Connect to Server` enter: `sftp://sftp.tudelft.nl/` and click `Connect`.
  - [ ] If requested, log in with your personal tudelft credentials.
  - [ ] Navigate to your data file (staff-bulk is the N-drive and staff-groups is the M-drive).
  - [ ] Copy the data set to the desktop, for example to `Documents/satay/datasets/`.
  - [ ] Close `Files`.
- [ ] Open the `Terminal` (located in the left bar) and enter `cd ~/Documents/satay/software/`.
- [ ] To start the workflow, enter `bash satay.sh`. This should open a window to select a fastq file.
- [ ] At the bottom right, select the right extension to be able to find your file, navigate to the location of your data set and select the fastq file and press `OK`. (For the test fastq, navigate to `Documents/satay/datasets/singleendtestfolder/SRR062634.fastq.gz`).
- [ ] In the next window (options window) enter the parameters and options. Click `Open adapters file` to change the sequences that need to be trimmed (if any) (this file should always be in fasta format, see [How to use](#how-to-use) section) and save and close the adapters file. Press `OK` to start processing.
- [ ] If the option `quality check interrupt` was selected:
  - [ ] Wait for the processing to finish the quality checking of the raw data.
  - [ ] The workflow will ask if you want to continue. Enter `y` to continue with the parameters you have set and enter `n` to stop the workflow to check the quality report of the raw data.
  - [ ] Navigate to the folder where your data file is stored and open the .html file (in the `fastqc_out` folder) and the check the quality report.
  - [ ] When finished, restart the workflow by typing `bash satay.sh` at `~/Documents/satay/software/`.
  - [ ] The workflow should skip the file selection window and immediately start with the processing options window where the parameters you have previously entered are set by default. Change the parameters if needed and press `OK` to continue processing.
  - [ ] The workflow should skip over the raw quality checking and continue with the processing.
- [ ] When processing is finished, navigate to folder where your dataset is located and check if all expected files are present (depending on what options yout have set, but at least the `align_out` folder and a log file should be present, see the [Input, Output](#input-output) section).
- [ ] Copy the results to the N-drive following the same procedure as described in step 2. (Be aware to prevent the sftp connection to automatically disconnect, see the note at the beginning of this section).
- [ ] After the processing is finished, most information is stored in the [bed](#bed) and [wig](#wig) files. But there can be some artifacts present in these files which can be removed using [strip_redundant_insertions.py](#strip_redundant_insertionspy). This is only necessary for specific downstream analysis tools like the [genome browser](#genome-browser).

#### How does it work

When the pipeline is started after the files and options are selected, it starts with checking all the paths for the software tools and data files.
The paths to the software tools and some other files can be adjusted by opening the bash script and change the paths in the `DEFINE PATHS` section at the beginning of the script.
Also it checks if the options that were given by the user don't conflict (e.g. when two data files were selected, but the option for processing the data as single-end reads in which case the second data file is ignored).
Next it defines names for all output files and folders.

After this the actual processing starts.
Here the full processing workflow is explained, but some parts may be skipped depending on which options are set by the user.

The first step is quality checking of the raw data in the fastq files using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
This creates a quality report in the `fastqc` folder and will take a few minutes depending on the size of the datafile.
After this first quality report the script will pause and ask the user to continue (enter `y` or `n`).
If the program is stopped, the quality report can be checked after which the program can be restarted (again with the command `bash satay.sh`).
All the previous settings will be remembered and do not have to be entered again.
But some options (e.g. trimming and alignment) can be altered if this is prefered by the user after checking the quality report.  
Next the trimming is started (either with [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/) or [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)).
If the option for paired-end data is selected, automatically these options are also selected for the trimming (e.g. `interleaved=t` for BBDuk or `PE` for Trimmomatic). Thus this should not be set manually by the user.
This step also reads the adapter sequences file to search for unwanted sequences in the reads.
The trimming can take more than half an hour to an hour (potentially even more for large datasets) and the results are stored in the `trimm_out` folder.
After the trimming, another quality report is created for the trimmed data which stored in the `fastqc` folder.  
Next up is the alignment by [BWA MEM](http://bio-bwa.sourceforge.net/bwa.shtml).
Just like with the trimming, arguments for indicating whether paired-end data is used (`-p`) is automatically set and therefore the user should not set these parameters manually.
The alignment can require more than one hour to complete and the results are stored in the `align_out` folder.  
After the alignment a quality report is created using [SAMTools](http://www.htslib.org/) which checks how many reads were aligned and if there are any things to be aware of (e.g. read pairs that were not mapped as proper pairs).
The results are stored as `_flagstatreport.txt` in the `align_out` folder.
The sam file is then converted to a bam file and it is checked whether the general format and structure of the bam file is correct.
Both tasks are done using SAMTools.
Next the bam file is sorted and an index file is created using [Sambamba](https://lomereiter.github.io/sambamba/).
All the steps so far after the the alignment should not take more than a few minutes and all the results are stored in the `align_out` folder.  
As a final step the python script [transposonmapping_satay.py](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/satay_processing/python_transposonmapping/transposonmapping_satay.py) is started (see below for more details).
This creates all the files that store the information about the insertions and reads in the `align_out` folder (i.e. bed, wig, pergene.txt, peressential.txt, pergene_insertions.txt and peressential_insertions.txt).
This python script can take well over an hour to complete.  
After the transposon mapping the sam file is deleted and a log file is created that stores information about which file(s) has been processed, a time stamp and all the options and adapter sequences set by the user.
This log file is stored together with the input fastq file.

The python script transposonmapping.py consist of a single python function called `transposonmapper`.
This takes a sorted bam file (which should have a bam index file present in the same folder) and potentially some additional information files.
It starts with checking if all input files are present.
Next if reads the bam file using the `pysam` function and gets some basic information about the chromosomes like the names and lengths and some statistics about the mapped reads in the bam file.  
Then a for-loop is started over all chromosomes and for each chromosome a for-loop is started over all reads.
For each read the samflag and the CIGAR string is taken to determine how the read was mapped (e.g. the orientation of the read).
By default the leftmost position is taken as insertion site, but if the read is reversely mapped, the length of the mapped part of the read needs to be added to the insertion position as the actual insertion happened at the rightmost position of the read.
After this correction the start position and the flags of each read are stored in arrays which are used later for creating the output files.  
Lists are created for all genes and all essential genes including their start and end positions in the genome.  
Next the chromosomes are all concatenated in one large array, so that each chromosome does not start anymore with basepair position 1, but rather continues counting from the previous chromosome (e.g. chrI starts at basepair 1 and ends at basepair 230218, then chrII starts at basepair 230219 instead of basepair 1).
This makes the next a bit easier where for each gene the insertions and reads are found that are mapped within the genes.
The last part of code consists of creating the different output files.
First the bed file is generated where the reads are saved using the equation (read_count*20)+100.
Then the pergene and peressential text files are created where the number of insertions and read count are stored together with the standard names of the genes.
The same is done for the pergene_insertions and peressential_insertions text files where all individual insertions and corresponding reads are stored which can be used later for determining the distribution of insertions within genes.
Finally the wiggle file is created where first the insertions that were mapped to the same location but with a different orientation are taken as a single insertion rather than two unique insertions.
The corresponding reads are added up and this information is stored as a wig file.

#### Notes

- [ ] The alignment software arguments all require to start with a dash (`-`). For some reason this sometimes gives errors and to prevent this make sure there is a space before the first argument in the `Enter alignment settings` line.
- [ ] If the workflow unexpectedly skips over the file selection window and immediately shows the options window, it might be that a previous processing run was not correctly finished. Press `cancel` and restart the workflow. It should now start with the file selection window. (This most likely has to do with a cachefile, that is created after the `quality check interrupt`, that was not deleted. This can be deleted manually as well by removing `processing_workflow_cache.txt` at the same location where the workflow is stored).
- [ ] When using the commandline for entering the arguments, use `Paired-end` and `Single-end` with first capital letters.
- [ ] For all software in the pipeline, whether the data is paired-end or single-end is automatically set, so do not set options like `interleaved=t` (for BBDuk), `PE` (for Trimmomatic)
or `-p` (for BWA) in the options window of the workflow.

## Software - analysis

### python scripts

Order to which to do the processing

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

For those who are not familiar with Linux, here is a brief summary of how to work with the Linux desktop and the basic tools that are necessary for using the workflow in Linux.

When using the Linux desktop, log in with the credentials found at `N:\tnw\BN\LL\Shared\LinuxMachines\LinuxDesktop_LoginCredentials.txt`.
For remote access, download [Teamviewer](https://www.teamviewer.com/nl/) on your computer and connect with the Linux desktop using the credentials found in the file mentioned above.

The workflow can be used via the commandline.
For this, open the Terminal located on the left bar in Linux.
This allows you to navigate through the different files and folders that are located on the desktop.
These files and folders can also be accessed using Files (also located on the left bar) which gives a more Windows like experience with clicking on pictograms, but this is missing some useful features from the commandline.
The Terminal by default opens in the home directory (indicated either as `~/` or `/home/laanlab/`).
For navigating through different directories (folders) use the change directory command `cd` (e.g. `cd ~/Documents/satay/` would bring you to the folder called satay which is located within the Documents folder).
Everything needed for the satay analysis is located in `~/Documents/satay`, which contains the following folders (among others):

- `software`; containing all the software packages and the bash script `satay.sh`.
- `datasets`; here the datasets are stored, including some test datasets.
- `reference_sequences`; containing the reference genome sequences.

While navigating through the folders, see what is inside a folder using `ls` (e.g. just typing `ls` shows the contents of the current directory and `ls ~/Documents/satay` would show the contents of the satay folder regardless in which folder you currently are).
There are several extensions of this command, like `ls -a` which shows also all hidden folders (starting with `.`) and `ls -l` shows some basic information about each file and folder.

To show the contents of a file use `less` or `head`.
The `head` command shows the first 10 lines of the file (e.g. `head ~/Documents/satay/software/satay.sh` shows the first 10 lines of the satay.sh script) and `less` also shows the contents of a file, but allows you to scroll through the file using the arrow keys.
The advantage of these commands are that they only load the lines that are being shown, so this is a perfect tool for checking large files without using a lot of memory (what would normally happen when opening a file completely).

To normally open a file with the default program for that filetype, use the command `xdg-open` (e.g. `xdg-open ~/Documents/satay/software/satay.sh` opens a text editor with satay.sh).

To access the webdrives, open the Files program and go to `Other Locations` (in the left menu bar in Files).
On the bottom of the window it should say `Connect to server`.
Enter here the following adress: `sftp://sftp.tudelft.nl/ `
Press `Connect` and, if requested, log in with your TUDelft credentials.
This should show a number of folders ,among of which is `staff-bulk` which is the N-drive and `staff-groups` which is the M-drive.
Here you can up-/download any files to the desktop, but after about 10 minutes of inactivity, the connection to the drives is automatically broken.
This means that all the up-/downloads will abort.
To prevent this, the switch back and forth to some files within the drives every few minutes to reset the timer.

As mentioned before, the workflow can be accessed by going to the following path: `~/Documents/satay/software/`
Start the workflow by stating the tool to which to run the script with (bash) and then the name of the script (satay.sh), thus entering `bash satay.sh`.
Any commands you want to put in come right after this statement, e.g. `bash satay.sh -h` to open the help text.
Datafiles can be stored at the following location: `~/Documents/satay/datasets/`.
To keep organized create a new folder for each dataset using the command `mkdir foldername` (where foldername is to be replaced with any name you want).
Prevent spaces in namegiving, but rather use underscores.
(When using the Files program, right mouse click and select `New Folder` at the location where the new folder has to be located).

Note that this desktop has no backup system, thus make sure to always put any files directly on the drives and remove the files you do not need from the desktop to prevent filling up the memory.

Once in a while, the desktop may require an update.
Sometimes the desktop will ask for this and you can follow the instructions of the pop up window.
To update manually enter the following commands in the commanline:

- `sudo apt-get update`; This checks for the latest versions of different system tools and files. This command only update the list of things to update, but doesn't do any downloading or isntalling itself. This requires a user password (required by the `sudo` command) which is the same password for logging in. (Note that when typing a password in the commandline, it doesn't show any signs that you're typing, but just enter the password and press enter).
  
- `sudo apt-get upgrade`; This does the actual downloading and installation of the updates that are found by the previous command. You may ask for a confirmation before it starts the installation and may restart the desktop.

For a more thorough tutorial of the linux commandline, see for example this course from [datacarpentry.org/shell-genomics](https://datacarpentry.org/shell-genomics/).
To have commandline like tool on windows (e.g. to access large datafiles on windows using `less` and `head` commands), try [git bash](https://gitforwindows.org/).
For Mac users, the terminal that is by default installed in MacOS works very similar as the one in Linux, so there is no real need for downloading any special tools.

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
