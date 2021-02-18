# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:06:48 2021

@author: gregoryvanbeek
"""

import os, sys
import numpy as np
from scipy import stats

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


read_gene_a = read_gene_df_a.loc[:, 'Nreadsperinsrt'].to_numpy()
read_gene_b = read_gene_df_b.loc[:, 'Nreadsperinsrt'].to_numpy()


#%% Define fold change and create new dataframe.
volcano_df = read_gene_df_a[['gene_names']]
#IN CASE OF WARNING: https://towardsdatascience.com/explaining-the-settingwithcopywarning-in-pandas-ebc19d799d25 . Add .copy() after .loc[...]
volcano_df['fold_change'] = read_gene_a - read_gene_b


#%% Determine t-statistic

#sample size
N = len(read_gene_df_a)

#variance
var_a = read_gene_a.var(ddof=1)
var_b = read_gene_b.var(ddof=1)

#standard_deviation
s = np.sqrt((var_a + var_b)/2)

#t-statistic
t = (read_gene_a.mean() - read_gene_b.mean()) / (s*np.sqrt(2/N))

#degrees of freedom
df = 2*N - 2

#p-value
p = 1 - stats.t.cdf(t, df=df)