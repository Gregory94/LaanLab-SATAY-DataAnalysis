# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 08:25:34 2020

@author: gregoryvanbeek

THIS SCRIPT ANALYSE THE INSERTION LOCATIONS PER GENE AS STORED IN _PERGENE_INSERTIONS.TXT AND PERESSENTIAL_INSERTIONS.TXT
"""

#%%
import os, sys
import copy
import pandas as pd
import seaborn as sns
import numpy as np


dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from gene_names import gene_aliases


#%%
def tninserts_analysis():
    '''
    IMPORTANT VARIABLES IN THIS FUNCTION:

    gene_position_dict: position of all genes
    gene_inserts_dict: insertion locations of transposons for all genes
    gene_reads_dict: number of reads for all tn inserts for all genes

    essential_position_dict: only for essential genes.
    essential_inserts_dict: only for essential genes.
    essential_reads_dict: only for essential genes.

    nonessential_position_dict: only for nonessential genes.
    nonessential_inserts_dict: only for nonessential genes.
    nonessential_reads_dict: only for nonessential genes.
    '''
#%% READ FILE AND PUT ALL VALUES IN DICTIONARIES
    filepath = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out"
    filename = "ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"
    datafile = os.path.join(filepath, filename)

    with open(datafile) as f:
        lines = f.readlines()


    gene_position_dict = {}
    gene_inserts_dict = {}
    gene_reads_dict = {}
    
    gene_inserts_distance_dict = {} #distance between subsequent inserts
    gene_inserts_trunc_dict = {} #inserts in the gene where 10% of the edges is truncated (so, only the center part of the gene is considered).
    for line in lines[1:]:
        line_split = line.strip('\n').split('\t')
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


        ins_list = []
        for ins in geneinserts_list: #GET INSERTIONS THAT ARE MORE THAN 10% OF LENGTH GENE AWAY FROM THE EDGES OF THE GENE.
            l = gene_end - gene_start
            if (gene_start + 0.1*l) < ins < (gene_end - 0.1*l):
                ins_list.append(ins)
        gene_inserts_trunc_dict[genename] = ins_list


        if not len(geneinserts_list) < 2:
            d = []
            for i in range(1,len(geneinserts_list)): #DISTANCES BETWEEN SUBSEQUENT INSERTS
                d.append(geneinserts_list[i] - geneinserts_list[i-1])
            gene_inserts_distance_dict[genename] = d
        elif len(geneinserts_list) == 1:
            gene_inserts_distance_dict[genename] = np.nan#[0] #only one insert
        else:
            gene_inserts_distance_dict[genename] = np.nan#[-1] #no insert


        genereads_str = line_split[5].strip('[]')
        if not genereads_str == '':
            genereads_list = [int(read) for read in genereads_str.split(',')]
        else:
            genereads_list = []
        gene_reads_dict[genename] = genereads_list


        if len(geneinserts_list) != len(genereads_list):
            print('Gene %s has different number of reads compared with the number of inserts' % genename )



    del (datafile, lines, line, line_split, genename, gene_chrom, gene_start, gene_end, geneinserts_str, geneinserts_list, genereads_str, genereads_list, i, d, ins, ins_list)
    #remains: gene_inserts_dict, gene_position_dict, gene_reads_dict


#%% GET ANNOTATED ESSENTIAL GENES
    essentialsfile = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_AllEssentialGenes_List.txt"

    with open(essentialsfile) as f:
        lines = f.readlines()


    aliases_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Yeast_Protein_Names.txt")[0]


    essential_position_dict = {}
    essential_inserts_dict = {}
    essential_reads_dict = {}

    nonessential_position_dict = copy.deepcopy(gene_position_dict)
    nonessential_inserts_dict = copy.deepcopy(gene_inserts_dict)
    nonessential_reads_dict = copy.deepcopy(gene_reads_dict)


    for line in lines[1:]:
        essential = line.strip('\n')


        essentiality = 'nonessential'


        if essential in gene_position_dict:
            essentiality = 'essential'
            alias = essential
        else:
            for alias in aliases_dict.get(essential):
                if alias in gene_position_dict:
                    essentiality = 'essential'
                    break


        if essentiality == 'essential':
            essential_position_dict[alias] = gene_position_dict.get(alias)
            essential_inserts_dict[alias] = gene_inserts_dict.get(alias)
            essential_reads_dict[alias] = gene_reads_dict.get(alias)
            
            del nonessential_position_dict[alias]
            del nonessential_inserts_dict[alias]
            del nonessential_reads_dict[alias]


    del (essentialsfile, lines, line, aliases_dict, essential, essentiality, alias)
    #remain: essential_position_dict, essential_inserts_dict, essential_reads_dict, nonessential_position_dict, nonessential_inserts_dict, nonessential_reads_dict



#%% CREATE DATAFRAME FOR ALL GENES
    genename_list = []
    essentiality_list = []
    N_inserts_list = []
    N_inserts_trunc_list = []
    distance_max_inserts_list = []
    N_reads_list = []
    for gene in gene_position_dict:
        genename_list.append(gene) #GENENAME LIST
        
        if gene in essential_position_dict: #ESSENTIALITY_LIST
            essentiality_list.append(True)
        elif gene in nonessential_position_dict:
            essentiality_list.append(False)
        else:
            print('WARNING: %s not found.' % gene)

        N_inserts_list.append(len(gene_inserts_dict.get(gene))) #N_INSERTS_LIST (NUMBER OF INSERTIONS)

        N_inserts_trunc_list.append(len(gene_inserts_trunc_dict.get(gene)))# / (gene_position_dict.get(gene)[2] - gene_position_dict.get(gene)[1])) #N_INSERTS_CENTER_LIST (NUMBER OF INSERTIONS IN THE GENE WHERE THE 10% OF THE GENE LENGTH IS TRUNCATED NORMALIZED TO GENE LENGTH)

        distance_max_inserts_list.append(np.max(gene_inserts_distance_dict.get(gene)) / (gene_position_dict.get(gene)[2] - gene_position_dict.get(gene)[1])) #DISTANCE_MAX_INSERTS_LIST (LARGEST DISTANCE BETWEEN SUBSEQUENT INSERTIONS NORMALIZED TO GENE LENGTH)
        
        N_reads_list.append(sum(gene_reads_dict.get(gene))) #N_READS_LIST



    allgenes = {'Gene_Name': genename_list,
                'Essentiality': essentiality_list,
                'Number_Insertions_Full_Gene': N_inserts_list,
                'Number_Insertions_Truncated_Gene': N_inserts_trunc_list,
                'Max_Insertion_Distance': distance_max_inserts_list,
                'Number_Reads': N_reads_list}


    df = pd.DataFrame(allgenes, columns = [column_name for column_name in allgenes])


#    del (gene, genename_list, essentiality_list, N_inserts_list, N_inserts_center_list, distance_max_inserts_list, N_reads_list, allgenes)


#%%TEST GRAPH
#    sns.set(style="whitegrid")
#    vp = sns.violinplot(x='Essentiality',y='Number_Insertions_Truncated_Gene',data=df, cut=0)
#    boxp = sns.boxplot(x='Essentiality',y='Number_Insertions_Truncated_Gene',data=df)

    bp = sns.barplot(x='Essentiality',y='Max_Insertion_Distance', data=df)

#    vp = sns.violinplot(x='Essentiality',y='Number_Reads', data=df, cut=0)




#%%
if __name__ == '__main__':
    tninserts_analysis()









