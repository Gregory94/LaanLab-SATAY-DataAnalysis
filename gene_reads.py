# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(1,r'C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\python_modules')
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic, gene_position
from gene_names import gene_aliases

#%%

def gene_reads(gene_name=None,region=None,wig_file=None):
    
#%% GET START AND END POSITION OF GENE
    if gene_name != None:
        gene_pos_dict = gene_position(r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3")
        
        gene_name = gene_name.upper() #CAPITALIZE GENE NAME
        if gene_pos_dict.get(gene_name) == None: #CHECK IF GENE_NAME EXISTS IN GENE_POS_DICT. IF NOT, CHECK IF ANY OF THE ALIASES EXISTS
            gene_alias_dict = gene_aliases()[0]
            gene_alias_key = [k for k, v in gene_alias_dict.items() if gene_name in v]
            print('gene_alias_key ',gene_alias_key[0])
            if gene_pos_dict.get(gene_alias_key[0]) == None: #IF KEY DOES ALSO NOT EXISTS IN GENE_POS_DICT, CHECK IF MORE ALIASES EXISTS OF GENE_NAME
                gene_alias_list = gene_alias_dict.get(gene_alias_key[0])
                for gene_alias in gene_alias_list:
                    if gene_pos_dict.get(gene_alias) != None:
                        gene_pos = gene_pos_dict.get(gene_alias)
                        print('The alias ',gene_alias, ' is used for the gene ',gene_name)
            else:
                gene_pos = gene_pos_dict.get(gene_alias_key[0])
                print('The alias ',gene_alias_key[0], ' is used for the gene ',gene_name)
                
        else:
            gene_pos = gene_pos_dict.get(gene_name)
        
    elif region != None:
        gene_pos = region
        

    gene_chr = gene_pos[0]
    gene_start = int(gene_pos[1])
    gene_end = int(gene_pos[2])
    if gene_name != None:
        print(gene_name, ' starts at basepair ',gene_start, ' and ends at basepair ',gene_end, ' in chromosome',gene_chr)
    else:
        print('Selected region starts at basepair ',gene_start, ' and ends at basepair ',gene_end, ' in chromosome',gene_chr)

#%% READ THE WIG FILE
    with open(wig_file) as f:
        lines = f.readlines()

#%% CREATE LIST OF CROMOSOME NUMBERS IN ROMAN NUMERALS
    arabic_to_roman_dict, roman_to_arabic_dict = chromosomename_roman_to_arabic()    
    chromosomenames_list = []
    for roman in roman_to_arabic_dict:
        chromosomenames_list.append(roman)

#%% GET THE NAMES OF THE CHROMOSOMES AS USED IN THE WIG FILE.
    chrom_names_dict = {} #STORES THE NAMES OF THE CHROMOSOMES
    chrom_lines_dict = {} #STORES WHERE THE CHROMOSOMES START IN THE WIG FILE
    chrom_names_counter = 0
    lines_counter = 0
    for line in lines:
        line.strip('\n')
        chrom_line = 'variableStep'
        line_split = line.split(' ')
        if line_split[0] == chrom_line:
            chromosome_name_wigfile = line_split[1].replace('chrom=chr','').strip('\n')
            chrom_names_dict[chromosomenames_list[chrom_names_counter]] = chromosome_name_wigfile
            chrom_lines_dict[chromosome_name_wigfile] = lines_counter
            print('Chromosome ',chromosomenames_list[chrom_names_counter], 'is ',chromosome_name_wigfile)
            
            chrom_names_counter += 1
        lines_counter += 1
    
#%% GET ALL LINES IN THE WIG FILE THAT INCLUDES THE WANTED CHROMOSOME
    chromosomenames_list_index = chromosomenames_list.index(gene_chr)
    chrom_range_start = chrom_lines_dict.get(chrom_names_dict.get(chromosomenames_list[chromosomenames_list_index]))+1
    if gene_chr == chromosomenames_list[-1]:
        chrom_range_end = len(lines)-1 #IF THE LAST CHROMOSOME IS NEEDED, GET ALL LINES TILL END OF FILE
    else:
        chrom_range_end = chrom_lines_dict.get(chrom_names_dict.get(chromosomenames_list[chromosomenames_list_index+1]))-1

#%% GET ALL READS WITHIN THE GENE
    insertion_list = []
    read_list = []
    for line in lines[chrom_range_start:chrom_range_end]:
        line_list = line.strip(' \n').split()
        if gene_start <= int(line_list[0]) <= gene_end:
            insertion_list.append(int(line_list[0]))
            read_list.append(int(line_list[1]))

#%% MAKE LIST OF ALL LOCATIONS IN THE GENE WITH THE NUMBER OF READS IN EACH LOCATION
    gene_length = gene_end-gene_start
    print('Length of region of interest is ',gene_length)
    roi_list = list(range(gene_start,gene_end+1))
    reads_roi_list = list(np.zeros(gene_length+1))
    
    read_index = 0
    for position in insertion_list:
        roi_index = roi_list.index(position)
        reads_roi_list[roi_index] = float(read_list[read_index])
        read_index += 1
    
#%% MAKE BAR PLOT FOR READS IN CHROMOSOME
    if gene_name != None:
        print('Plotting reads for gene ', gene_name, '...')
    else:
        print('Plotting reads in range ', gene_start, '..', gene_end, 'in chromosome ', gene_chr, '...')
    fig,ax = plt.subplots()
    ax.bar(roi_list,reads_roi_list,width=10)
    ax.set_axisbelow(True)
    ax.grid(True)
    if gene_name != None:
        ax.set_title(gene_name)
    ax.set_xlabel('Basepair position in chromosome '+ gene_chr, fontsize=12)        
    ax.set_ylabel('Read count', fontsize=12)
    ax.set_xlim(gene_start,gene_end)
    plt.show()

#%%
if __name__ == '__main__':
#    gene_reads(region=['IV',46271,48031],wig_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.wig")
    gene_reads(gene_name='snu71',wig_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.wig")