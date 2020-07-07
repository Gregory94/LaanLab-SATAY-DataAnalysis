# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import os, sys
import timeit

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import gene_position
from gene_names import gene_aliases
from essential_genes_names import list_known_essentials

#%% Using primarily for-loops

#timer_start = timeit.default_timer()
#
## GET POSITION GENES
#gff_path = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
#gene_pos_dict = gene_position(gff_path) #contains all genes, essential and nonessential
#
## GET ALL ANNOTATED ESSENTIAL GENES
#essential_path1 = os.path.join(file_dirname,'Data_Files','Cervisiae_EssentialGenes_List_1.txt')
#essential_path2 = os.path.join(file_dirname,'Data_Files','Cervisiae_EssentialGenes_List_2.txt')
#known_essential_gene_list = list_known_essentials([essential_path1, essential_path2])
#
## GET ALIASES OF ALL GENES
#names_path = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')
#aliases_designation_dict = gene_aliases(names_path)[0]
#
## FOR ALL ANNOTATED ESSENTIAL GENES, DETERMINE THEIR ALIASES AND GET THE POSITION FROM GENE_POS_DICT
#essential_pos_dict = {}
#for gene in gene_pos_dict:
#    if gene in known_essential_gene_list:
#        essential_pos_dict[gene] = gene_pos_dict.get(gene)
#    else:
#        gene_aliases_list = []
#        for key, val in aliases_designation_dict.items():
#            if gene == key or gene in val: #if gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
#                gene_aliases_list.append(key)
#                for aliases in aliases_designation_dict.get(key):
#                    gene_aliases_list.append(aliases)
#    
#        for gene_alias in gene_aliases_list:
#                if gene_alias in known_essential_gene_list:
#                    essential_pos_dict[gene_alias] = gene_pos_dict.get(gene)
#                    break
#
#timer_end = timeit.default_timer()
#print('Running time is %.2f seconds' %float(timer_end-timer_start))


#%% Using map() function

timer_start = timeit.default_timer()

def search_aliases_dict(gene):
    gene_aliases_list = []
    if gene == key or gene in val: #if gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
        gene_aliases_list.append(key)
        for aliases in aliases_designation_dict.get(key):
            gene_aliases_list.append(aliases)
    return(gene_aliases_list)




# GET POSITION GENES
gff_path = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
gene_pos_dict = gene_position(gff_path) #contains all genes, essential and nonessential

# GET ALL ANNOTATED ESSENTIAL GENES
essential_path1 = os.path.join(file_dirname,'Data_Files','Cervisiae_EssentialGenes_List_1.txt')
essential_path2 = os.path.join(file_dirname,'Data_Files','Cervisiae_EssentialGenes_List_2.txt')
known_essential_gene_list = list_known_essentials([essential_path1, essential_path2])

# GET ALIASES OF ALL GENES
names_path = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')
aliases_designation_dict = gene_aliases(names_path)[0]

# FOR ALL ANNOTATED ESSENTIAL GENES, DETERMINE THEIR ALIASES AND GET THE POSITION FROM GENE_POS_DICT
essential_pos_dict = {}
for gene in gene_pos_dict:
    if gene in known_essential_gene_list:
        essential_pos_dict[gene] = gene_pos_dict.get(gene)
    else:
        gene_aliases_list = map(search_aliases_dict, gene)
    
        for gene_alias in gene_aliases_list:
                if gene_alias in known_essential_gene_list:
                    essential_pos_dict[gene_alias] = gene_pos_dict.get(gene)
                    break

timer_end = timeit.default_timer()
print('Running time is %.2f seconds' %float(timer_end-timer_start))