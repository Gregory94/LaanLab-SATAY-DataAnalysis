# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:38:47 2021

@author: gregoryvanbeek
"""
import os

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