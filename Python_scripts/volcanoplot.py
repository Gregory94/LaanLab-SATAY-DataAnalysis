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


#%% define file paths and names. Two samples called a and b.
datapath_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing"
filenames_a = ["WT-a_pergene.txt", "WT-b_pergene.txt"]
datapath_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing"
filenames_b = ["dnrp1-1-a_pergene.txt", "dnrp1-1-b_pergene.txt"]


#%% Check files
datafiles_list_a = []
datafiles_list_b = []
for files in filenames_a:
    datafile = os.path.join(datapath_a, files)
    assert os.path.isfile(datafile), 'File not found at: %s' % datafile
    datafiles_list_a.append(datafile)
for files in filenames_b:
    datafile = os.path.join(datapath_b, files)
    assert os.path.isfile(datafile), 'File not found at: %s' % datafile
    datafiles_list_b.append(datafile)

del (files, datafile, datapath_a, datapath_b, filenames_a, filenames_b)



#%% Get information from files and store in pandas dataframe
#ADD COLUMN IN READPERINSRT_DF WITH THE STANDARD DEVIATION
counter = 0
for datafile_a in datafiles_list_a:
    read_gene_a = dataframe_from_pergenefile(datafile_a)
    if counter == 0:
        readsperinsrt_df = read_gene_a[['gene_names']] #initialize new dataframe with gene_names
        readsperinsrt_df['mean_Nreadsperinsrt_a'] = read_gene_a[['Nreadsperinsrt']] #add column with Nreadsperinsrt
        Nreadspersinsrt_a_array = read_gene_a[['Nreadsperinsrt']].to_numpy() #create numpy array to store raw data
    else:
        readsperinsrt_df['mean_Nreadsperinsrt_a'] = readsperinsrt_df['mean_Nreadsperinsrt_a'] + read_gene_a['Nreadsperinsrt'] #add values of Nreadsperinsrt to existing Nreadsperinsrt column
        Nreadspersinsrt_a_array = np.append(Nreadspersinsrt_a_array, read_gene_a[['Nreadsperinsrt']].to_numpy(), axis=1) #append raw data
    counter += 1
readsperinsrt_df['mean_Nreadsperinsrt_a'] = readsperinsrt_df['mean_Nreadsperinsrt_a'].div(counter) #determine mean Nreadsperinsrt
N_a = counter #save number of samples



counter = 0
for datafile_b in datafiles_list_b:
    read_gene_b = dataframe_from_pergenefile(datafile_b)
    if counter == 0:
        readsperinsrt_df['mean_Nreadsperinsrt_b'] = read_gene_b[['Nreadsperinsrt']]
        Nreadspersinsrt_b_array = read_gene_b[['Nreadsperinsrt']].to_numpy()
    else:
        readsperinsrt_df['mean_Nreadsperinsrt_b'] = readsperinsrt_df['mean_Nreadsperinsrt_b'] + read_gene_b['Nreadsperinsrt']
        Nreadspersinsrt_b_array = np.append(Nreadspersinsrt_b_array, read_gene_b[['Nreadsperinsrt']].to_numpy(), axis=1)
    counter += 1
readsperinsrt_df['mean_Nreadsperinsrt_b'] = readsperinsrt_df['mean_Nreadsperinsrt_b'].div(counter)
N_b = counter



del (datafile_a, datafile_b, counter, read_gene_a, read_gene_b)


#%% Get information from files and store in pandas dataframe BACKUP
#IN CASE OF WARNING: https://towardsdatascience.com/explaining-the-settingwithcopywarning-in-pandas-ebc19d799d25 . Add .copy() after .loc[...]
# counter = 0
# for datafile_a in datafiles_list_a:
#     read_gene_a = dataframe_from_pergenefile(datafile_a)
#     if counter == 0:
#         readsperinsrt_df = read_gene_a[['gene_names']]
#         readsperinsrt_mean_array = read_gene_a[['Nreadsperinsrt']].to_numpy()
#     readsperinsrt_df['mean_Nreadsperinsrt_a'+str(counter)] = read_gene_a[['Nreadsperinsrt']]
#     readsperinsrt_mean_array = np.add(read_gene_a[['Nreadsperinsrt']].to_numpy(), readsperinsrt_mean_array)
#     if counter == len(datafiles_list_a):
#         readsperinsrt_mean_array = readsperinsrt_mean_array / np.array(counter)
#     counter += 1


# counter = 0
# for datafile_b in datafiles_list_b:
#     read_gene_b = dataframe_from_pergenefile(datafile_b)
#     readsperinsrt_df['mean_Nreadsperinsrt_b'+str(counter)] = read_gene_b[['Nreadsperinsrt']]
#     counter += 1


# del (datafile_a, datafile_b, read_gene_a, read_gene_b)


#%% Define fold change and create new dataframe.
#IN CASE OF WARNING: https://towardsdatascience.com/explaining-the-settingwithcopywarning-in-pandas-ebc19d799d25 . Add .copy() after .loc[...]
# volcano_df = readsperinsrt_df[['gene_names']]
# volcano_df['fold_change'] = readsperinsrt_df - readsperinsrt_df












