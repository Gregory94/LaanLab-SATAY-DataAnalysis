# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import sys
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
sys.path.insert(1,r'C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\python_modules')
from chromosome_and_gene_positions import chromosomename_roman_to_arabic, gene_position
from gene_names import gene_aliases

#%%

def gene_reads(gene_name=None,region=None,bed_file=None):
    '''This script makes a profile plot for the number of reads per tranposon for a specific genomic region.
    Input is a region and the .bed file from the output of the Matlab code from the Kornmann-lab.
    The region can be defined either as a gene name (e.g. 'bem1') or as a list consisting of three elements where the first element is the chromosome name, the start and end position respectively (e.g. ['I',1,4000]).
    If a gene name is input, the script searches in a .gff file (downloaded from yeastgenome.org).
    The output is a bar plot where the number of reads divided by the number of transposons.
    '''
#%% USED FILES
    gff_file = r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    gene_information_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt'

#%% GET START AND END POSITION OF GENE
    if gene_name != None:
        gene_pos_dict = gene_position(gff_file)
        
        gene_name = gene_name.upper() #CAPITALIZE GENE NAME
        if gene_pos_dict.get(gene_name) == None: #CHECK IF GENE_NAME EXISTS IN GENE_POS_DICT. IF NOT, CHECK IF ANY OF THE ALIASES EXISTS
            gene_alias_dict = gene_aliases(gene_information_file)[0]
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
    with open(bed_file) as f:
        lines = f.readlines()

#%% CREATE LIST OF CROMOSOME NUMBERS IN ROMAN NUMERALS
    arabic_to_roman_dict, roman_to_arabic_dict = chromosomename_roman_to_arabic()    
    chromosome_romannames_list = []
    for roman in roman_to_arabic_dict:
        chromosome_romannames_list.append(roman)

#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
    chrom_names_dict = {}
    chrom_start_line_dict = {}
    chrom_end_line_dict = {}
    chrom_name_in_bed = ''
    chr_counter = 0
    line_counter = 0
    stop_loop = False
    while stop_loop is False:
        line = lines[line_counter]
        chrom_name_current = line.split(' ')[0].replace('chr','')
        if not chrom_name_current.startswith('track') and not chrom_name_current.startswith('M'): #SKIP HEADER AND MITOCHRONDRIAL CHROMOSOMES
            if chrom_name_current != chrom_name_in_bed:
                chrom_names_dict[chromosome_romannames_list[chr_counter]] = chrom_name_current
                chrom_name_in_bed = chrom_name_current
                print('Chromosome ',chromosome_romannames_list[chr_counter], 'is ',chrom_name_current)
                
                chrom_start_line_dict[chromosome_romannames_list[chr_counter]] = line_counter #GET START INDEX IN THE BED FILE OF THE CURENT CHROMOSOME
                if chr_counter != 0:
                    chrom_end_line_dict[chromosome_romannames_list[chr_counter-1]] = line_counter-1 #GET THE END INDEX IN THE BED OF THE PREVIOUS CHROMOSOME (SKIP FOR THE FIRST CHROMOSOME)

                chr_counter += 1

        elif chrom_name_current.startswith('M'):
            chrom_end_line_dict[chromosome_romannames_list[-1]] = line_counter-1 #GET THE END INDEX IN THE BED FILE FOR THE FINAL CHROMOSOME
            stop_loop = True
                
        line_counter += 1

#%% GET ALL READS WITHIN THE GENE
    insertion_list = []
    read_list = []
    for line in lines[chrom_start_line_dict.get(gene_chr):chrom_end_line_dict.get(gene_chr)]:
        line_list = line.strip('\n').split()
        if gene_start <= int(line_list[1]) <= gene_end:
            insertion_list.append(int(line_list[1]))
            
            read_value = (int(line_list[4])-100)/20
            read_list.append(read_value)

#%% ACCOUNT FOR DOUBLE INSERTIONS
#see chromosome I, bp 3891
    unique_insertion_list = []
    duplicate_insertion_list = []
    for ins in insertion_list: #FIND ALL DUPLICATED TRANSPOSON INSERTION SITES
        if ins not in unique_insertion_list:
            unique_insertion_list.append(ins)
        else:
            duplicate_insertion_list.append(ins)
    duplicate_insertion_list = np.unique(duplicate_insertion_list)
    
    duplicate_index_list = []
    for dup in duplicate_insertion_list:
        insertion_arr = np.asarray(insertion_list)
        duplicate_index_list.append(np.where(insertion_arr == dup))
    
    if len(duplicate_index_list) > 0:
        number_of_duplicates_list = [1]*len(insertion_list) #MAKE LIST OF ONES WITH SAME LENGTH AS INSERTION_LIST FOR STORING NUMBER OF DUPLICATES
        delete_index = []
        for ind_arr in duplicate_index_list:#LOOP OVER ALL INDICES OF DUPLICATES
            ind_list = ind_arr[0]
            ind_list_max = max(ind_list)#GET THE LAST INDEX OF THE DUPLICATES
            print('Mulitple transposons found at ',ind_list)
            for ind in ind_list:
                if not ind == ind_list_max:
                    read_list[ind_list_max] += read_list[ind]#ADD UP THE READS TO THE LAST DUPLICATE
                    number_of_duplicates_list[ind_list_max] = len(ind_list) #UPDATE NUMBER OF DUPLICATES
                    delete_index.append(ind)
        
        #REVERSE LOOP OVER LIST FOR DELETING
        for del_ind in reversed(delete_index):
            del read_list[del_ind] #DELETES THE INDEX WHICH IS NOW ADDED UP TO THE LAST INDEX
            del insertion_list[del_ind] #DELETES THE SAME INDICES IN THE INSERTION LISTS.
            del number_of_duplicates_list[del_ind]
        
        readspertransposon_list = [x/y for x, y in zip(read_list, number_of_duplicates_list)] #DIVIDE THE NUMBER OF READS BY THE NUMBER OF TRANSPOSONS
    else:
        readspertransposon_list = read_list

#%% MAKE LIST OF ALL LOCATIONS IN THE GENE WITH THE NUMBER OF READS IN EACH LOCATION
    gene_length = gene_end-gene_start
    print('Length of region of interest is ',gene_length)
    roi_list = list(range(gene_start,gene_end+1))
    reads_roi_list = list(np.zeros(gene_length+1))
#    color_bars_roi_list = list(np.zeros(gene_length+1))
    
    read_index = 0
    for position in insertion_list:
        roi_index = roi_list.index(position)
        reads_roi_list[roi_index] = float(readspertransposon_list[read_index])
#        color_bars_roi_list[roi_index] = float(number_of_duplicates_list[read_index])
        read_index += 1
    
#%% MAKE BAR PLOT FOR READS IN CHROMOSOME

#    clist = [(1, "green"), (max(color_bars_roi_list), "red")]
#    rvb = mcolors.LinearSegmentedColormap.from_list("", clist)

    if gene_name != None:
        print('Plotting reads for gene ', gene_name, '...')
    else:
        print('Plotting reads in range ', gene_start, '..', gene_end, 'in chromosome ', gene_chr, '...')
    fig,ax = plt.subplots()
    ax.bar(roi_list,reads_roi_list,width=15,color='k')#, color=rvb(roi_list/len(roi_list)))
    ax.set_axisbelow(True)
    ax.grid(True)
    if gene_name != None:
        ax.set_title(gene_name, fontweight='bold')
    elif region == ['IV',46271,48031]:
        ax.set_title('HO-locus', fontweight='bold')
    ax.set_xlabel('Basepair position in chromosome '+ gene_chr)        
    ax.set_ylabel('Read/tn')
    ax.set_xlim(gene_start,gene_end)
    plt.show()

#%%
if __name__ == '__main__':
    gene_reads(region=['IV',46271,48031],bed_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed")
#    gene_reads(gene_name='bem3',bed_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed")