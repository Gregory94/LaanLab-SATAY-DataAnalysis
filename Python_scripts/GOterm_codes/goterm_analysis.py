# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:57:07 2020

@author: gregoryvanbeek

THIS SCRIPT CREATES AN OVERVIEW OF ALL GO-TERMS WITH THE CORRESPONDING GENES.
THE GO-TERMS WITH THE GENES IS DOWNLOADED FROM http://sgd-archive.yeastgenome.org/curation/literature/
THE LIST OF GENES IS USED FROM S_Cerevisiae_protein_oln_name_full_genome.txt CREATED USING gene_names.py
"""

#%%
import pandas as pd


#%% LOAD FILES
gene_file = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\S_Cerevisiae_protein_oln_name_full_genome.txt"
go_file = r"C:\Users\gregoryvanbeek\Documents\GeneOntology\go_slim_mapping_SGD.tab"

#%% LOAD ALL GENES FROM GENE_FILE
with open(gene_file, 'r') as f:
    lines = f.readlines()
    
gene_list = []
for line in lines[1:]:
    gene_list.append(line.strip('\n'))


del (lines, line, f)


#%% READ GO_FILE
with open(go_file,'r') as f:
    lines = f.readlines()

go_dict = {}
for line in lines:
    l = line.strip('\n').split('\t')
    if not l[0] in go_dict:
        go_dict[l[0]] = [l[1:]]
    else:
        go_dict[l[0]].append(l[1:])

del (lines, line, l, f)


#%% CREATE DICTIONARIES WITH GO TERMS AND GENE NAMES

for gene in gene_list:
    
