# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import os, sys
import numpy as np

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname))
from chromosome_and_gene_positions import gene_position
from gene_names import gene_aliases
from chromosome_names_in_files import chromosome_name_bedfile

#%%

def hit_free_region(gene_name='None',region=None,bed_file=None):
    '''This script makes a profile plot for the number of reads per tranposon for a specific genomic region.
    Input is a region and the .bed file from the output of the Matlab code from the Kornmann-lab.
    The region can be defined either as a gene name (e.g. 'bem1') or as a list consisting of three elements where the first element is the chromosome name, the start and end position respectively (e.g. ['I',1,4000]).
    If a gene name is input, the script searches in a .gff file (downloaded from yeastgenome.org).
    The output is a bar plot where the number of reads divided by the number of transposons.
    '''
#%% USED FILES
    datafile_dirname = os.path.join(file_dirname,'..')
    
    gff_file = os.path.join(datafile_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    gene_information_file = os.path.join(datafile_dirname,'Data_Files','Yeast_Protein_Names.txt')

#%% GET START AND END POSITION OF GENE
    if gene_name.upper() == 'HOLOCUS' or gene_name == 'HO-LOCUS':
        gene_pos = ['IV',46271,48031]
        gene_name = 'HOlocus'

    elif gene_name != 'None':
        gene_pos_dict = gene_position(gff_file) #GET POSITION INFORMATION OF ALL GENES
        
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


#%% READ THE BED FILE
    with open(bed_file) as f:
        lines = f.readlines()

#%% GET POSITION FOR THE CHROMOSOMES IN THE BED FILE

    chrom_start_line_dict, chrom_end_line_dict= chromosome_name_bedfile(lines)[1:3]

#%% GET ALL READS WITHIN THE GENE
    insertion_list = []
    read_list = []
    for line in lines[chrom_start_line_dict.get(gene_chr):chrom_end_line_dict.get(gene_chr)]:
        line_list = line.strip('\n').split()
        if gene_start <= int(line_list[1]) <= gene_end:
            insertion_list.append(int(line_list[1]))

            read_value = (int(line_list[4])-100)/20 #the matlab script by benoit takes the number of reads*20+100. This line makes this undone
            read_list.append(read_value)

#%% ACCOUNT FOR DOUBLE INSERTIONS FOR PLOTTING
#see for example chromosome I, bp 3891
    unique_insertion_list = []
    duplicate_insertion_list = []
    for ins in insertion_list: #FIND THE CHROMOSOME POSITION OF ALL DUPLICATED TRANSPOSON INSERTION SITES
        if ins not in unique_insertion_list:
            unique_insertion_list.append(ins)
        else:
            duplicate_insertion_list.append(ins)
    duplicate_insertion_list = np.unique(duplicate_insertion_list) #ACCOUNT FOR THE SITUATION WHERE THERE ARE MORE THAN TWO INSERTIONS AT THE SAME LOCATION

    duplicate_index_list = []
    for dup in duplicate_insertion_list:
        insertion_arr = np.asarray(insertion_list)
        duplicate_index_list.append(np.where(insertion_arr == dup)) #GET ALL INDICES OF THE LIST OF TRANSPOSON INSERTIONS WHERE THE DUPLICATES ARE PRESENT. EACH INSERTION LOCATION GETS ITS OWN NUMPY ARRAY WITHIN THIS LIST

    if len(duplicate_index_list) > 0:
        number_of_duplicates_list = [1]*len(insertion_list) #MAKE LIST OF ONES WITH SAME LENGTH AS INSERTION_LIST FOR STORING NUMBER OF DUPLICATES
        delete_index = []
        for ind_arr in duplicate_index_list: #LOOP OVER ALL INDICES OF DUPLICATES
            ind_list = ind_arr[0]
            ind_list_max = max(ind_list) #GET THE LAST INDEX OF THE DUPLICATES
#            print('Mulitple transposons found at ',ind_list)
            for ind in ind_list:
                if not ind == ind_list_max:
                    read_list[ind_list_max] += read_list[ind] #ADD UP THE READS TO THE LAST DUPLICATE
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
    insertion_roi_list = list(range(gene_start,gene_end+1))
    reads_roi_list = list(np.zeros(gene_length+1))
    
    read_index = 0
    for position in insertion_list:
        roi_index = insertion_roi_list.index(position)
        reads_roi_list[roi_index] = float(readspertransposon_list[read_index])
        read_index += 1

#%% CALCULATE SOME STATISTICAL VALUES FOR THE SELECTED REGION
#insertion_roi_list := list of all potential insertion sites in the region
#reads_roi_list := number of reads in the selected region.

    if insertion_list != []:
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(insertion_list[:-1], insertion_list[1:])]
        bp_between_tn_insertions.insert(0,insertion_list[0] - gene_start) #ADD START OF GENE (bp=0)
        bp_between_tn_insertions.append(gene_end - insertion_list[-1]) #ADD END OF GENE (bp=INDEX LAST TN - GENE LENGTH)

        max_empty_region = max(bp_between_tn_insertions)
        
    else:
        max_empty_region = gene_length
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(insertion_list[:-1], insertion_list[1:])]

#%%
    print('insertion_list: ', insertion_list)
    print('read_list: ', read_list)
    print('Basepairs between subsequent insertions: ', bp_between_tn_insertions) #FIRST AND LAST DIFFERENCE IS NUMBER OF BASEPAIRS BETWEEN FIRST AND LAST INSERTION AND THE BEGINNING AND END OF THE REGION
    print('max_empty_region: ', max_empty_region)
    
    return(insertion_list, read_list, max_empty_region, bp_between_tn_insertions)

#%%
if __name__ == '__main__':
#    hit_free_region(region=['XI',1000,2000],bed_file=r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\satay_analysis_testdata\Output_Processing\Cerevisiae_WT2_Michel2017_trimmed1.bam.bed")
    hit_free_region(gene_name='cdc24',bed_file=r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\satay_analysis_testdata\Output_Processing\Cerevisiae_WT2_Michel2017_trimmed1.bam.bed")
