# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:38:47 2021

@author: gregoryvanbeek
"""
import os
from gene_names import gene_aliases

gene_alias_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Data_Files\Yeast_Protein_Names.txt")[0]


filepath_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\leila_dataset_wt_processing\WT_merged-techrep-a_techrep-b_processing\WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"
filepath_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\WT_pergene.txt"

assert os.path.isfile(filepath_a), "Filepath_a not found."
assert os.path.isfile(filepath_b), "Filepath_b not found."

genes_a = []
with open(filepath_a, 'r') as a:
    lines = a.readlines()[1:]
for line in lines:    
    if '\t' in line:
        genes_a.append(line.strip('\n').split('\t')[0])
    else:
        genes_a.append(line.strip('\n').split(' ')[0])


genes_b = []
with open(filepath_b, 'r') as b:
    lines = b.readlines()[1:]
for line in lines:    
    if '\t' in line:
        genes_b.append(line.strip('\n').split('\t')[0])
    else:
        genes_b.append(line.strip('\n').split(' ')[0])

del (a, b, line, lines)



genes_in_a_not_in_b = []
genes_in_b_not_in_a = []
for gene in genes_a:
    if not gene in genes_b:
        genes_in_a_not_in_b.append(gene)
for gene in genes_b:
    if not gene in genes_a:
        genes_in_b_not_in_a.append(gene)

del (gene)


#check for aliases
for index, gene in enumerate(genes_in_a_not_in_b):
    alias_list = []
    for keys, nums in gene_alias_dict.items():
        if keys == gene or gene in nums:
            alias_list = gene_alias_dict.get(keys)
            alias_list.append(keys)
            break

    if not alias_list == []:
        for alias in alias_list:
            if alias in genes_in_b_not_in_a:
                genes_in_b_not_in_a.remove(alias)
                break


# del (keys, nums, gene, alias_list, index, alias)
