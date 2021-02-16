# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:56:11 2021

@author: gregoryvanbeek
"""
#%% initialization
import os
import pandas as pd

filepath = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing\WT-a_pergene.txt"

assert os.path.isfile(filepath), 'File not found at: %s' % filepath


#%% read file
with open(filepath) as f:
    lines = f.readlines()[1:] #skip header

genenames_list = [None]*len(lines)
tnpergene_list = [None]*len(lines)
readpergene_list = [None]*len(lines)  

line_counter = 0
for line in lines:
    l = line.strip('\n').split(' ')

    genenames_list[line_counter] = l[0]
    tnpergene_list[line_counter] = l[1]
    readpergene_list[line_counter] = l[2]

    line_counter += 1

del (lines, line, l, line_counter, filepath)


#%% determine essential genes



#%% create dataframe

all_features = {"gene_names": genenames_list,
                "tn_per_gene": tnpergene_list,
                "read_per_gene": readpergene_list}

genes_df = pd.DataFrame(all_features, columns = [column_name for column_name in all_features])

del (all_features)
