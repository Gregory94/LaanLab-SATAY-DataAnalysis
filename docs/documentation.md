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

This documentation gives a complete overview for the processing of the data from SAturated Transposon Analysis in Yeast (SATAY).
A short introduction to SATAY will be given and details how to perform the processing is discussed in more detail.

For this a pipeline is created using Bash and Python.
The workflow and the python codes can be found at [github.com/Gregory94/LaanLab-SATAY-DataAnalysis](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/satay_processing).
More information about satay analysis and experimental protocols can be found at the [satayusers website from the Kornmann lab](https://sites.google.com/site/satayusers/ "satayusers website") or, for more questions, visit the [satayusers forum](https://groups.google.com/g/satayusers "satayusers forum").

## Introduction

## File types

### fastq

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
