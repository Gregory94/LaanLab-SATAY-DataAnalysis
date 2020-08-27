# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:02:44 2020

@author: gregoryvanbeek

This scripts takes a user defined genomic region (i.e. chromosome number, region or gene) and determines the number of transposon insertions in the genomic features (i.e. genes, nc-DNA etc.).
This can be used to determine the number of transposon insertions outside the genes to use this for normalization of the number of transposons in the genes.
"""

#%%
import os, sys
import numpy as np
import pandas as pd


file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic
from chromosome_names_in_files import chromosome_name_wigfile


#%% USER INPUT AND FILES

region = 'II'

wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533148_trimmed.sorted.bam.wig"

pergene_insertions_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"

gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')


#%% DETERMINE INPUTTED REGION

if type(region) == str:
    if region in chromosomename_roman_to_arabic()[1]:
        chrom = region
    else:
        print("WARNING: Specified chromosome not found. Enter chromosome as a roman numeral between I and XVI")

elif type(region) == list:
    chrom = region[0]
    region_start = region[1]
    region_end = region[2]


del (region)

#%% READ WIG FILE FOR GETTING ALL TN INSERTIONS

with open(wig_file, 'r') as f:
    lines = f.readlines()

chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(lines)[1:]

insrt_in_chrom_list = []
reads_in_chrom_list = []
for l in lines[chrom_start_line_dict.get(chrom):chrom_end_line_dict.get(chrom)]:
    insrt_in_chrom_list.append(int(l.strip('\n').split(' ')[0]))
    reads_in_chrom_list.append(int(l.strip('\n').split(' ')[1]))


del (lines, l, f, chrom_start_line_dict, chrom_end_line_dict)

#%% READ PERGENE_INSERTIONS FILE FOR LOCATION AND INSERTIONS FOR EACH GENE.

with open(pergene_insertions_file) as f:
    lines = f.readlines()


gene_position_dict = {}
gene_inserts_dict = {}
gene_reads_dict = {}
for line in lines[1:]:
    line_split = line.strip('\n').split('\t')


    if line_split[1] == chrom:
        genename = line_split[0]
        gene_chrom = line_split[1]
        gene_start = int(line_split[2])
        gene_end = int(line_split[3])
    
        gene_position_dict[genename] = [gene_chrom, gene_start, gene_end]
    
    
        geneinserts_str = line_split[4].strip('[]')
        if not geneinserts_str == '':
            geneinserts_list = [int(ins) for ins in geneinserts_str.split(',')]
        else:
            geneinserts_list = []
    
        gene_inserts_dict[genename] = geneinserts_list
    
    
        genereads_str = line_split[5].strip('[]')
        if not genereads_str == '':
            genereads_list = [int(read) for read in genereads_str.split(',')]
        else:
            genereads_list = []
    
        gene_reads_dict[genename] = genereads_list
    
    
        if len(geneinserts_list) != len(genereads_list):
            print('WARNING: %s has different number of reads compared with the number of inserts' % genename )


del (f, lines, line, genename, gene_chrom, gene_start, gene_end, geneinserts_list, geneinserts_str, genereads_str, genereads_list)

#%% DETERMINE THE LOCATION GENOMIC FEATURES IN THE CURRENT CHROMOSOME AND STORE THIS IN A DICTIONARY

start_chr = chromosome_position(gff_file)[1].get(chrom)
end_chr = chromosome_position(gff_file)[2].get(chrom)

dna_dict = {}
for bp in range(start_chr, end_chr + 1):
    dna_dict[bp] = 'noncoding'


for gene in gene_position_dict: #IGNORE GENES THAT END WITH -A AS THESE ARE DUBIOUS (?) 
    for bp in range(gene_position_dict.get(gene)[1]+start_chr, gene_position_dict.get(gene)[2]+start_chr+1):
        dna_dict[bp] = gene
        
#coding_dna_bplocations_dict = {}
#N_bp_in_genes = 0 #needed for preallocating coding_dna_bplocations_array with the right length
#for gene in gene_position_dict: #GENE_POSITION_DICT ALSO INCLUDES DUBIOUS OPEN READING FRAMES.
#    coding_dna_bplocations_dict[gene] = list(range(gene_position_dict.get(gene)[1]+start_chr, gene_position_dict.get(gene)[2]+start_chr+1))
#    N_bp_in_genes += (gene_position_dict.get(gene)[2] - gene_position_dict.get(gene)[1]) + 1
#
#coding_dna_bplocations_array = np.zeros(N_bp_in_genes) #contains all unique basepairs at location of gene
#counter = 0
#for gene in coding_dna_bplocations_dict:
#    for i in coding_dna_bplocations_dict.get(gene):
#        coding_dna_bplocations_array[counter] = int(i)
#        counter += 1
#
#coding_dna_bplocations_array.sort()
#coding_dna_bplocations_array = np.unique(coding_dna_bplocations_array).astype(int)


#ALL LOCATIONS THAT DO NOT OCCUR IN CODING_DNA_BPLOCATIONS_LIST.
#noncoding_dna_bplocations_array = np.array(range(start_chr, end_chr+1), dtype=int)
#del_bp_ind_list = [0] * ((end_chr - start_chr) - len(coding_dna_bplocations_array) + 3)
#i = 0
#j = 0
#for v in noncoding_dna_bplocations_array:
#    if v in coding_dna_bplocations_array:
#        del_bp_ind_list[j] = i
#        j += 1
#    i += 1
#
#noncoding_dna_bplocations_array = np.delete(noncoding_dna_bplocations_array, del_bp_ind_list)


#noncoding_dna_bplocations_array = np.zeros((end_chr - start_chr) - len(coding_dna_bplocations_array) + 2, dtype=int)
#counter = 1
#for j in range(start_chr, end_chr+1):
#    if not int(j) in coding_dna_bplocations_array: #SEARCHING IS SLOW
#        noncoding_dna_bplocations_array[counter] = int(j)
#        counter += 1



del (gene, bp)

#%% DETERMINE THE NUMBER OF TRANSPOSONS PER BP FOR EACH FEATURE AND STORE RESULTS IN DATAFRAME

dna_df = pd.DataFrame(list(dna_dict.items()), columns=["BP", "Feature"])


reads_loc_list = [0] * len(dna_dict) # CONTAINS ALL READS JUST LIKE READS_IN_CHROM_LIST, BUT THIS LIST HAS THE SAME LENGTH AS THE NUMBER OF BP IN THE CHROMOSOME WHERE THE LOCATIONS WITH NO READS ARE FILLED WITH ZEROS
i = 0
for ins in insrt_in_chrom_list:
    reads_loc_list[ins] = reads_in_chrom_list[i]
    i += 1
dna_df["Reads"] = reads_loc_list


del (i, ins)

#%% CREATE DATAFRAME FOR EACH FEATURE (E.G. NONCODING DNA, GENE, ETC.) IN THE CHROMOSOME AND DETERMINE THE NUMBER OF INSERTIONS AND READS PER FEATURE.
feature_name_list = []
f_previous = dna_dict.get(start_chr)
N_reads = 0
N_reads_list = []
N_insrt = 0
N_insrt_list = []
N_bp = 0
N_bp_list = []
i = 0
for bp in dna_dict:
    f_current = dna_dict.get(bp)
    if f_current == f_previous:
        N_bp += 1
        N_reads += reads_loc_list[i]
        if not reads_loc_list[i] == 0:
            N_insrt += 1
    elif f_current != f_previous or (i+start_chr) == end_chr:
        feature_name_list.append(f_previous)
        N_reads_list.append(N_reads)
        N_insrt_list.append(N_insrt)
        N_bp_list.append(N_bp)
        N_reads = 0
        N_insrt = 0
        N_bp = 0
        f_previous = f_current
    i += 1

N_reads_per_bp_list = []
N_insrt_per_bp_list = []
for i in range(len(N_reads_list)):
    N_reads_per_bp_list.append(N_reads_list[i]/N_bp_list[i])
    N_insrt_per_bp_list.append(N_insrt_list[i]/N_bp_list[i])

all_features = {'Feature': feature_name_list,
                'Nreads':N_reads_list,
                'Ninsertions':N_insrt_list,
                'Nbasepairs':N_bp_list,
                'Nreadspertn':N_reads_per_bp_list,
                'Ninsertionspertn':N_insrt_per_bp_list}

dna_df2 = pd.DataFrame(all_features, columns = [column_name for column_name in all_features])




del (feature_name_list, f_previous, f_current, N_reads, N_reads_list, N_insrt, N_insrt_list, N_bp, N_bp_list, bp, i, N_reads_per_bp_list, N_insrt_per_bp_list, all_features)


#%% PLOTTING

#SEE TRANSPOSONREAD_PROFILE_PLOT.PY; READ_PROFILE FUNCTION