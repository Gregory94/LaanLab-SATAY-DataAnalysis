#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a tool developed for analysing transposon insertions for experiments using SAturated Transposon Analysis in Yeast (SATAY).
For more information about this code and the project, see github.com/Gregory94/LaanLab-SATAY-DataAnalysis

This code is based on the Matlab code created by the Kornmann lab and is available at: sites.google.com/site/satayusers/
It contains one function called transposonmapper().

__Author__ = Gregory van Beek. LaanLab, department of Bionanoscience, Delft University of Technology
__version__ = 1.0
__Date last update__ = 2020-07-24
"""

import os, sys
import warnings
import numpy as np
import pysam

dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from chromosome_and_gene_positions import chromosomename_roman_to_arabic, gene_position
from gene_names import gene_aliases
from essential_genes_names import list_known_essentials


def transposonmapper(bamfile=None, gfffile=None, essentialfiles=None, genenamesfile=None):
    '''
    This function is created for analysis of SATAY data using the species Saccharomyces Cerevisiae.
    It outputs the following files that store information regarding the location of all insertions:
        - .bed-file: Include all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Browser Extensible Data (bed) format.
                    A distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored using the equation #reads*20+100 (e.g. 2 reads is stored as 140).
        - .wig-file: Include all individual basepair locations of the whole genome where at least one transposon has been mapped and the number of insertions for each locations (the number of reads) according to the Wiggle (wig) format.
                    In this file no distinction is made between reads that had a different reading orientation during sequencing. The number of reads are stored as the absolute count.
        - _pergene.txt-file: Include all genes (currently 6600) with the total number of insertions and number of reads within the genomic region of the gene.
    The output files are saved at the location of the input file using the same name as the input file, but with the corresponding extension.
    
    The function assumes that the reads are already aligned to a reference genome.
    The input data should be a .bam-file and the location where the .bam-file is stored should also contain an index file (.bai-file, which for example can be created using sambamba).
    It takes the following inputs:
        - bamfile [required]: Path to the bamfile. This location should also contain the index file (.bai).
        - gfffile [optional]: Path to a .gff-file including all gene information (e.g. downloaded from SGD).
        - essentialfiles [optional]: List of paths with text files that includes all essential gene names.
        - genenamesfile [optional]: Path to text file that includes aliases for all genes.
    When the arguments for the optional files are not given, the files are used that are stored at the following location:
        "path_current_pythonscript/../data_files"
    The function uses the pysam package for handling bam files (see pysam.readthedocs.io/en/latest/index.html)
    '''

#%% LOADING BAM FILE
    if bamfile is None:
        path = os.path.join('/home', 'gregoryvanbeek', 'Documents', 'data_processing')
        filename = 'E-MTAB-4885.WT1.bam'
        bamfile = os.path.join(path,filename)
    
    if os.path.isfile(bamfile):
        print('Running: ', bamfile)
    else:
        raise ValueError('Bam file not found.')



#%% LOADING ADDITIONAL FILES
    files_path = os.path.join(dirname,'..','data_files')
    
    #LOADING GFF-FILE
    if gfffile is None:
        gfffile = os.path.join(files_path,'Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    if not os.path.isfile(gfffile):
        raise ValueError('Path to GFF-file does not exist.')
    
    #LOADING TEXT FILES WITH ESSENTIAL GENES
    if essentialfiles is None:
        essentialfiles = [os.path.join(files_path,'Cervisiae_EssentialGenes_List_1.txt'), os.path.join(files_path,'Cervisiae_EssentialGenes_List_2.txt')]
    for essential_file in essentialfiles:
        if not os.path.isfile(essential_file):
            raise ValueError('Following path does not exist: ' + essential_file)
    del essential_file
    
    #LOADING TEXT FILE WITH GENE NAME ALIASES
    if genenamesfile is None:
        genenamesfile = os.path.join(files_path,'Yeast_Protein_Names.txt')
    if not os.path.isfile(genenamesfile):
        raise ValueError('Following path does not exist: ' + genenamesfile)



#%% READ BAM FILE
    bam = pysam.AlignmentFile(bamfile, 'rb') #open bam formatted file for reading



#%% GET NAMES OF ALL CHROMOSOMES AS STORED IN THE BAM FILE
    ref_tid_dict = {} # 'I' | 0, 'II' | 1, ...
    ref_name_list = [] # 'I', 'II', ...
    for i in range(bam.nreferences): #if bam.nreferences does not work, use range(17) #16 chromosomes and the mitochondrial chromosome
        ref_name = bam.get_reference_name(i)
        ref_tid_dict[ref_name] = bam.get_tid(ref_name)
        ref_name_list.append(ref_name)

    del (ref_name, i)


#%% GET SEQUENCE LENGTHS OF ALL CHROMOSOMES
    chr_length_dict = {} # 'I' | 230218, 'II' | 813184, ...
    chr_summedlength_dict = {} # 'I' | 0, 'II' | 230218, ...
    ref_summedlength = 0
    for key in ref_tid_dict:
        ref_length = bam.get_reference_length(key)
        chr_length_dict[key] = ref_length
        chr_summedlength_dict[key] = ref_summedlength
        ref_summedlength += ref_length

    del (key, ref_length, ref_summedlength)



#%% GET NUMBER OF MAPPED, UNMAPPED AND TOTAL AMOUNT OF READS PER CHROMOSOME
    # total_reads = bam.mapped
    stats = bam.get_index_statistics()
    chr_mappedreads_dict = {} # 'I' | [mapped, unmapped, total reads]
    for stat in stats:
        chr_mappedreads_dict[stat[0]] = [stat[1], stat[2], stat[3]]
        if stat[2] != 0:
            warnings.warn('Unmapped reads found in chromosome ' + stat[0])
    
    del (stat, stats)



#%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION





#%%










#%%
if __name__ == '__main__':
    transposonmapper()



