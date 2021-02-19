# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:06:48 2021

@author: gregoryvanbeek

What we are doing when we do a t test is measuring the number of standard
deviations our measured mean is from the baseline mean, while taking into
account that the standard deviation of the mean can change as we get more data
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



#%% Get number of reads per insertion all datasets
for count, datafile_a in enumerate(datafiles_list_a):
    read_gene_a = dataframe_from_pergenefile(datafile_a)
    if count == 0:
        readsperinsrt_df = read_gene_a[['gene_names']] #initialize new dataframe with gene_names
        Nreadsperinsrt_a_array = read_gene_a[['Nreadsperinsrt']].to_numpy() #create numpy array to store raw data
    else:
        Nreadsperinsrt_a_array = np.append(Nreadsperinsrt_a_array, read_gene_a[['Nreadsperinsrt']].to_numpy(), axis=1) #append raw data
N_a = count+1 #save number of samples



for count, datafile_b in enumerate(datafiles_list_b):
    read_gene_b = dataframe_from_pergenefile(datafile_b)
    if count == 0:
        Nreadsperinsrt_b_array = read_gene_b[['Nreadsperinsrt']].to_numpy()
    else:
        Nreadsperinsrt_b_array = np.append(Nreadsperinsrt_b_array, read_gene_b[['Nreadsperinsrt']].to_numpy(), axis=1)
N_b = count+1


if not N_a == N_b:
    print("WARNING: Length of dataset a is NOT the same as the length of dataset b.")
N = count+1

del (datafile_a, datafile_b, count, read_gene_a, read_gene_b, N_a, N_b)




#%% Determine mean and p value of reads per insertion for each gene
df = 2*N - 2 #degrees of freedom

mean_Nreadsperinsrt_a_list = [np.nan]*len(Nreadsperinsrt_a_array)
mean_Nreadsperinsrt_b_list = [np.nan]*len(Nreadsperinsrt_b_array)
var_Nreadsperinsrt_a_list = [np.nan]*len(Nreadsperinsrt_a_array)
var_Nreadsperinsrt_b_list = [np.nan]*len(Nreadsperinsrt_b_array)

for ii in range(len(Nreadsperinsrt_a_array)): #determine mean and variance
    mean_Nreadsperinsrt_a_list[ii] = Nreadsperinsrt_a_array[ii].mean()
    mean_Nreadsperinsrt_b_list[ii] = Nreadsperinsrt_b_array[ii].mean()

    #For unbiased max likelihood estimate we have to divide the var by N-1, and therefore the parameter ddof = 1
    var_Nreadsperinsrt_a_list[ii] = Nreadsperinsrt_a_array[ii].var(ddof=1)
    var_Nreadsperinsrt_b_list[ii] = Nreadsperinsrt_b_array[ii].var(ddof=1)



std_Nreadsperinsrt_list = [np.nan]*len(Nreadsperinsrt_a_array)
tstat_Nreadsperinsrt_list = [np.nan]*len(Nreadsperinsrt_a_array)
p_Nreadsperinsrt_list = [np.nan]*len(Nreadsperinsrt_a_array)

for ss in range(len(var_Nreadsperinsrt_a_list)): #determine standard deviation and t statistic
    std_Nreadsperinsrt_list[ss] = np.sqrt((var_Nreadsperinsrt_a_list[ss] + var_Nreadsperinsrt_b_list[ss])/2)
    if std_Nreadsperinsrt_list[ss] == 0.0:
        tstat_Nreadsperinsrt_list[ss] = 0
    else:
        tstat_Nreadsperinsrt_list[ss] = (mean_Nreadsperinsrt_a_list[ss] - mean_Nreadsperinsrt_b_list[ss]) / (std_Nreadsperinsrt_list[ss] * np.sqrt(2/N))
    p_Nreadsperinsrt_list[ss] = 1 - stats.t.cdf(tstat_Nreadsperinsrt_list[ss], df=df)

readsperinsrt_df['mean_Nreadsperinsrt_a'] = mean_Nreadsperinsrt_a_list
readsperinsrt_df['mean_Nreadsperinsrt_b'] = mean_Nreadsperinsrt_b_list
readsperinsrt_df['std_Nreadsperinsrt'] = std_Nreadsperinsrt_list
readsperinsrt_df['tstat_Nreadsperinsrt'] = tstat_Nreadsperinsrt_list
readsperinsrt_df['pval_Nreadsperinsrt'] = p_Nreadsperinsrt_list


del (mean_Nreadsperinsrt_a_list, mean_Nreadsperinsrt_b_list, var_Nreadsperinsrt_a_list, var_Nreadsperinsrt_b_list, ii, ss,
     std_Nreadsperinsrt_list, tstat_Nreadsperinsrt_list, p_Nreadsperinsrt_list)


#%% Determine fold change of mean_Nreadsperinsrt
fc_list = [np.nan]*len(readsperinsrt_df) #initialize list for storing fold changes

for count, avg in enumerate(readsperinsrt_df.itertuples()):
    if not avg.mean_Nreadsperinsrt_b == 0:
        fc_list[count] = np.log2(avg.mean_Nreadsperinsrt_a / avg.mean_Nreadsperinsrt_b) #DIVIDE DATASET A BY DATASET B
    else:
        fc_list[count] = 0

readsperinsrt_df['log2_fold_change'] = fc_list #add fc_list to dataframe

del (avg, count, fc_list)


#%% Volcanoplot









