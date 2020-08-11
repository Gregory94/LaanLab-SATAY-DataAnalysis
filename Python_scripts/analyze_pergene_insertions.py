# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 08:25:34 2020

@author: gregoryvanbeek

THIS SCRIPT ANALYSE THE INSERTION LOCATIONS PER GENE AS STORED IN _PERGENE_INSERTIONS.TXT AND PERESSENTIAL_INSERTIONS.TXT
"""

#%%
import os, sys
import copy


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
    filepath = r"C:\Users\gregoryvanbeek\Documents\testing_site\WT2_dataset_analysis_temp202008051429\new2"
    filename = "E-MTAB-4885.WT2.bam_pergene_insertions.txt"
    datafile = os.path.join(filepath, filename)

    with open(datafile) as f:
        lines = f.readlines()


    gene_position_dict = {}
    gene_inserts_dict = {}
    gene_reads_dict = {}
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


        genereads_str = line_split[5].strip('[]')
        if not genereads_str == '':
            genereads_list = [int(read) for read in genereads_str.split(',')]
        else:
            genereads_list = []
        gene_reads_dict[genename] = genereads_list


        if len(geneinserts_list) != len(genereads_list):
            print('Gene %s has different number of reads compared with the number of inserts' % genename )


    del (datafile, lines, line, line_split, genename, gene_chrom, gene_start, gene_end, geneinserts_str, geneinserts_list, genereads_str, genereads_list)
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


#%% APPLY STATISTICS ON THE NINE DICTIONARIES
    for gene in gene_position_dict:
        #number of tn inserts
        #largest gap between tn inserts
        #histogram of distances between tn inserts
        #Make criteria for predicting which genes are essential. Check with annotation.
        #...



#%%
if __name__ == '__main__':
    tninserts_analysis()









