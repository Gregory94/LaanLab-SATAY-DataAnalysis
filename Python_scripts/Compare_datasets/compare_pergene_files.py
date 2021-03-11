# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:38:47 2021

@author: gregoryvanbeek
"""
import os

filepath_a = r"C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\SATAY\data\15022021-sequencing-data-WT-dnrp1-SATAY-from-Oxford\a-b-pooled\WT_pergene.txt"
filepath_greg = r"X:\tnw\BN\LL\Shared\Gregory\datasets\dataset_leila\leila_dataset\leila_dataset_wt\leila_dataset_wt_processing\WT_merged-techrep-a_techrep-b_processing\align_out\WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene.txt"

assert os.path.isfile(filepath_a), "Filepath_a not found."
assert os.path.isfile(filepath_greg), "Filepath_b not found."

genes_a = []
with open(filepath_a, 'r') as a:
    lines = a.readlines()[1:]
for line in lines:    
    if '\t' in line:
        genes_a.append(line.strip('\n').split('\t')[0])
    else:
        genes_a.append(line.strip('\n').split(' ')[0])


genes_b = []
with open(filepath_greg, 'r') as b:
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
#%% Finding aliases genes between those datasets
from read_sgdfeatures import sgd_features
from gene_names import  gene_aliases

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sgd_features_file = os.path.join(file_dirname,'..','..','Data_Files','SGD_features.tab')

gene_information_file = os.path.join(file_dirname,'..','..','Data_Files','Yeast_Protein_Names.txt')

gene_alias_dict = gene_aliases(gene_information_file)[0] # it is only used systematic gene names
feature_orf_dict = sgd_features(sgd_features_file)[1]

same_genes_for_a=[]
unique_genes_in_b_not_in_a=[]
unique_genes_in_a_not_in_b=[]
k=0
for gene in genes_in_a_not_in_b: # for the case of their analyses they sometimes use systematic gene names 
    
    if gene.startswith('Y') and len(gene)>4: # looking for systematic gene annotation
        if gene_alias_dict[gene][0] == genes_in_b_not_in_a[k]:
            same_genes_for_a.append(genes_in_b_not_in_a[k])
        else :
            unique_genes_in_b_not_in_a.append(genes_in_b_not_in_a[k])
            unique_genes_in_a_not_in_b.append(genes_in_a_not_in_b[k])
    
    
    
    k=k+1
    
k=0
same_genes_for_b=[]
for gene in genes_in_b_not_in_a: # for the case of their analyses they sometimes use systematic gene names 
    
    if gene.startswith('Y') and len(gene)>4: # looking for systematic gene annotation
        if feature_orf_dict[gene][2] == genes_in_a_not_in_b[k] :
            same_genes_for_b.append(genes_in_a_not_in_b[k])
        else :
            unique_genes_in_b_not_in_a.append(genes_in_b_not_in_a[k])
            unique_genes_in_a_not_in_b.append(genes_in_a_not_in_b[k])
    
    
    
    k=k+1