# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:06:48 2021

@author: gregoryvanbeek
"""

import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from dataframe_from_pergene import dataframe_from_pergenefile


datafile_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing\WT-a_pergene.txt"
datafile_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing\dnrp1-1-a_pergene.txt"


assert os.path.isfile(datafile_a), 'File not found at: %s' % datafile_a
assert os.path.isfile(datafile_b), 'File not found at: %s' % datafile_b


read_gene_df_a = dataframe_from_pergenefile(datafile_a)
read_gene_df_b = dataframe_from_pergenefile(datafile_b)

del (datafile_a, datafile_b, file_dirname)


assert (len(read_gene_df_a) == len(read_gene_df_b)), "Lengths dataframe_a and dataframe_b are not the same. Please check data."


fc = read_gene_df_a[['gene_names']]
fc['fold_change'] = read_gene_df_a.Nreadsperinsrt - read_gene_df_b.Nreadsperinsrt

