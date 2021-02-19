# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:06:48 2021

@author: gregoryvanbeek

This script creates a volcanoplot to show the significance of fold change between two datasets.
The following steps are taken:
    - Load all datafiles from two samples, called a and b
    - For each gene in each sample, determine the mean and the unbiased variance
    - For each gene determine the standard deviation based on the variances in each sample.
    - For each gene determine the t-statistic
    - For each gene determine the negative log_10 of the p-value corresponding to the t-statistic
    - For each gene determine the log_2 fold change between the means of number of reads per insertion of both samples: log2((mean_a-mean_b)/mean_b)
    - plot log fold change vs negative log p value
This is based on this website:
    - https://towardsdatascience.com/inferential-statistics-series-t-test-using-numpy-2718f8f9bf2f
    - https://www.statisticshowto.com/independent-samples-t-test/

T-test is measuring the number of standard
deviations our measured mean is from the baseline mean, while taking into
account that the standard deviation of the mean can change as we get more data
"""

import os, sys
import numpy as np
from scipy import stats
# import seaborn as sns
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from dataframe_from_pergene import dataframe_from_pergenefile



#%% define file paths and names. Two samples called a and b.
datapath_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing"
filenames_a = ["WT-a_pergene.txt", "WT-b_pergene.txt"]
datapath_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing"
filenames_b = ["dnrp1-1-a_pergene.txt", "dnrp1-1-b_pergene.txt", "dnrp1-2-a_pergene.txt", "dnrp1-2-b_pergene.txt"]



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



#%% Extract information from datasets
variable = 'read_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
print('Plotting: %s' % variable)

#!!! NORMALIZE FOR TOTAL NUMBER OF INSERTIONS?

for count, datafile_a in enumerate(datafiles_list_a):
    read_gene_a = dataframe_from_pergenefile(datafile_a)
    if count == 0:
        readsperinsrt_df = read_gene_a[['gene_names']] #initialize new dataframe with gene_names
        Nreadsperinsrt_a_array = read_gene_a[[variable]].to_numpy() #create numpy array to store raw data
    else:
        Nreadsperinsrt_a_array = np.append(Nreadsperinsrt_a_array, read_gene_a[[variable]].to_numpy(), axis=1) #append raw data
N_a = count+1 #save number of samples



for count, datafile_b in enumerate(datafiles_list_b):
    read_gene_b = dataframe_from_pergenefile(datafile_b)
    if count == 0:
        Nreadsperinsrt_b_array = read_gene_b[[variable]].to_numpy()
    else:
        Nreadsperinsrt_b_array = np.append(Nreadsperinsrt_b_array, read_gene_b[[variable]].to_numpy(), axis=1)
N_b = count+1


# if not N_a == N_b:
    # print("WARNING: Length of dataset a is NOT the same as the length of dataset b.")
N = count+1 #!!!THIS N IS NEEDED IN NEXT SECTION, BUT DOES NOT WORK WHEN DATASETS DO NOT INCLUDE THE SAME NUMBER OF DATA FILES.

del (datafile_a, datafile_b, count, variable) #, read_gene_a, read_gene_b)




#%% Determine mean and p value of reads per insertion for each gene
df = N_a + N_b - 2 #degrees of freedom

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
        #!!! Check tstat. WHAT TO DO WITH THE N AT THE END?
        tstat_Nreadsperinsrt_list[ss] = (mean_Nreadsperinsrt_a_list[ss] - mean_Nreadsperinsrt_b_list[ss]) / (std_Nreadsperinsrt_list[ss] * np.sqrt(2/N))
    p_Nreadsperinsrt_list[ss] = -1*np.log10(stats.t.cdf(tstat_Nreadsperinsrt_list[ss], df=df))
    # p_Nreadsperinsrt_list[ss] = 1 - stats.t.cdf(tstat_Nreadsperinsrt_list[ss], df=df)

readsperinsrt_df['mean_Nreadsperinsrt_a'] = mean_Nreadsperinsrt_a_list
readsperinsrt_df['mean_Nreadsperinsrt_b'] = mean_Nreadsperinsrt_b_list
readsperinsrt_df['std_Nreadsperinsrt'] = std_Nreadsperinsrt_list
readsperinsrt_df['tstat_Nreadsperinsrt'] = tstat_Nreadsperinsrt_list
readsperinsrt_df['pval_Nreadsperinsrt'] = p_Nreadsperinsrt_list


del (mean_Nreadsperinsrt_a_list, mean_Nreadsperinsrt_b_list, var_Nreadsperinsrt_a_list, var_Nreadsperinsrt_b_list, ii, ss,
     std_Nreadsperinsrt_list, tstat_Nreadsperinsrt_list, p_Nreadsperinsrt_list, df, N)



#%% Determine fold change of mean_Nreadsperinsrt
fc_list = [np.nan]*len(readsperinsrt_df) #initialize list for storing fold changes

for count, avg in enumerate(readsperinsrt_df.itertuples()):
    if not avg.mean_Nreadsperinsrt_a == 0 and not avg.mean_Nreadsperinsrt_b == 0:
        fc_list[count] = np.log2((avg.mean_Nreadsperinsrt_b / avg.mean_Nreadsperinsrt_a)) #DIVIDE DATASET A BY DATASET B
        # fc_list[count] = (avg.mean_Nreadsperinsrt_a / avg.mean_Nreadsperinsrt_b)-1
    else:
        fc_list[count] = np.log2(max([avg.mean_Nreadsperinsrt_a, avg.mean_Nreadsperinsrt_b]))

readsperinsrt_df['log2_fold_change'] = fc_list #add fc_list to dataframe

del (avg, count, fc_list)


#print specific geness
# readsperinsrt_df[readsperinsrt_df['gene_names'].str.contains("NRP1")]



#%% Volcanoplot
plt.figure(figsize=(19.0,9.0))#(27.0,3))
grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])

# sns.scatterplot(data=readsperinsrt_df, x='log2_fold_change', y='pval_Nreadsperinsrt', alpha=0.6)
ax.scatter(x=readsperinsrt_df.log2_fold_change, y=readsperinsrt_df.pval_Nreadsperinsrt, alpha=0.4, marker='.', c='k')
ax.grid(True, which='major', axis='both', alpha=0.4)
ax.set_xlabel('Log2 FC')
ax.set_ylabel('-Log10 p-value')



#%%THIS IS AN ALTERNATIVE APPROACH FOR THE ABOVE CALCULATION. START FROM "Extract information from datasets"
#%% TEST INDEPENDENT T-TEST
## https://www.statisticshowto.com/independent-samples-t-test/

# test1=[541, 664]
# test2=[799,396,711,567]

# len_test1 = len(test1)
# len_test2 = len(test2)

# sum_test1 = sum(test1)
# sum_test2 = sum(test2)

# mean_test1 = np.mean(test1)
# mean_test2 = np.mean(test2)

# sum_sqrt_test1 = 0
# sum_sqrt_test2 = 0
# for val in test1:
#     sum_sqrt_test1 += val**2
# for val in test2:
#     sum_sqrt_test2 += val**2

# t1 = mean_test1 - mean_test2
# t2 = (sum_sqrt_test1 - (sum_test1**2 / len_test1)) + (sum_sqrt_test2 - (sum_test2**2 / len_test2))
# t3 = len_test1 + len_test2 - 2
# t4 = (1/len_test1) + (1/len_test2)
# t = t1 / np.sqrt((t2/t3)*t4)

# print(t)
# print(stats.ttest_ind(test1,test2))

#%% APPLY stats.ttest_ind(A,B)













