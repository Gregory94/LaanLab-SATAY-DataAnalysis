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
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from dataframe_from_pergene import dataframe_from_pergenefile



#%% define file paths and names. Two samples called a and b.
datapath_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing"
filenames_a = ["WT-a_pergene.txt", "WT-b_pergene.txt"]
datapath_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing"
filenames_b = ["dnrp1-1-a_pergene.txt", "dnrp1-1-b_pergene.txt", "dnrp1-2-a_pergene.txt", "dnrp1-2-b_pergene.txt"]


variable = 'tn_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

track_gene = "yal069w"# "CDC42" or set to "" to disable


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
print('Plotting: %s' % variable)

# norm_a = 0
# norm_b = 0
for count, datafile_a in enumerate(datafiles_list_a):
    tnread_gene_a = dataframe_from_pergenefile(datafile_a, verbose=False)
    if normalize == True:
        if variable == 'tn_per_gene':
            norm_a = sum(tnread_gene_a.tn_per_gene)#*10**-4
        elif variable == 'read_per_gene':
            norm_a = sum(tnread_gene_a.read_per_gene)#*10**-7
    # norm_a += sum(tnread_gene_a.tn_per_gene)

    if count == 0:
        volcano_df = tnread_gene_a[['gene_names']] #initialize new dataframe with gene_names
        if normalize == True:
            variable_a_array = np.divide(tnread_gene_a[[variable]].to_numpy(), norm_a) #create numpy array to store raw data
        else:
            variable_a_array = tnread_gene_a[[variable]].to_numpy()
    else:
        if normalize == True:
            variable_a_array = np.append(variable_a_array, np.divide(tnread_gene_a[[variable]].to_numpy(), norm_a), axis=1) #append raw data
        else:
            variable_a_array = np.append(variable_a_array, tnread_gene_a[[variable]].to_numpy(), axis=1)


for count, datafile_b in enumerate(datafiles_list_b):
    tnread_gene_b = dataframe_from_pergenefile(datafile_b, verbose=False)
    if normalize == True:
        if variable == 'tn_per_gene':
            norm_b = sum(tnread_gene_b.tn_per_gene)#*10**-4
        elif variable == 'read_per_gene':
            norm_b = sum(tnread_gene_b.read_per_gene)#*10**-7
    # norm_b += sum(tnread_gene_b.tn_per_gene)

    if count == 0:
        if normalize == True:
            variable_b_array = np.divide(tnread_gene_b[[variable]].to_numpy(), norm_b)
        else:
            variable_b_array = tnread_gene_b[[variable]].to_numpy()
    else:
        if normalize == True:
            variable_b_array = np.append(variable_b_array, np.divide(tnread_gene_b[[variable]].to_numpy(), norm_b), axis=1)
        else:
            variable_b_array = np.append(variable_b_array, tnread_gene_b[[variable]].to_numpy(), axis=1)


### printing specific genes
if not track_gene == "":
    track_gene = track_gene.upper()
    print("TRACKING VARIABLE %s" % track_gene)
    track_gene_index = tnread_gene_a.loc[tnread_gene_a['gene_names'] == track_gene].index[0]
    print("tnread_gene_a:", tnread_gene_a.loc[tnread_gene_a['gene_names'] == track_gene])
    print("tnread_gene_b:", tnread_gene_b.loc[tnread_gene_b['gene_names'] == track_gene])


del (datafile_a, datafile_b, count, tnread_gene_a, tnread_gene_b)


#%% APPLY stats.ttest_ind(A,B)
fc_list = [np.nan]*len(variable_a_array) #initialize list for storing fold changes
ttest_tval_list = [np.nan]*len(variable_a_array) #initialize list for storing t statistics
ttest_pval_list = [np.nan]*len(variable_a_array) #initialize list for storing p-values
signif_thres_list = [False]*len(variable_a_array) #initialize boolean list for indicating datapoints with p-value above threshold

for count, val in enumerate(variable_a_array):

    ttest_val = stats.ttest_ind(variable_a_array[count], variable_b_array[count]) #T-test
    ttest_tval_list[count] = ttest_val[0]
    if not ttest_val[1] == 0: #prevent p=0 to be inputted in log
        ttest_pval_list[count] = -1*np.log10(ttest_val[1])
    else:
        ttest_pval_list[count] = 0
    if ttest_pval_list[count] > -1*np.log10(significance_threshold):
        signif_thres_list[count] = True

    #Take mean of number of insertions per library
    if np.mean(variable_b_array[count]) == 0 and np.mean(variable_a_array[count]) == 0:
        fc_list[count] = 0
    elif np.mean(variable_b_array[count]) == 0 or np.mean(variable_a_array[count]) == 0:
        fc_list[count] = np.log2(max(np.mean(variable_a_array[count]), np.mean(variable_b_array[count])))
    else:
        fc_list[count] = np.log2(np.mean(variable_b_array[count]) / np.mean(variable_a_array[count]))

    #Take sum of number of insertions per library 
    # if sum(variable_a_array[count]) == 0 and sum(variable_b_array[count]) == 0:
    #     fc_list[count] = 0
    # elif sum(variable_a_array[count]) == 0 or sum(variable_b_array[count]) == 0:
    #     fc_list[count] = np.log2(max(sum(variable_a_array[count]) / norm_a, sum(variable_b_array[count]) / norm_b))
    # else:
    #     fc_list[count] = np.log2((sum(variable_a_array[count]) / norm_a) / (sum(variable_b_array[count]) / norm_b))


    ### printing specific genes
    if not track_gene == "" and count == track_gene_index:
        print("variable_a_array:", variable_a_array[count])
        print("variable_b_index:", variable_b_array[count])
        print("t-test value:", ttest_val)
        print("mean a:", np.mean(variable_a_array[count]))
        print("mean b:", np.mean(variable_b_array[count]))
        print("fold change:", fc_list[count])


volcano_df['fold_change'] = fc_list
volcano_df['t_statistic'] = ttest_tval_list
volcano_df['p_value'] = ttest_pval_list
volcano_df['significance'] = signif_thres_list

del(count, val, ttest_val, ttest_tval_list, ttest_pval_list, fc_list, signif_thres_list)
if normalize == True:
    del (norm_a, norm_b)


#%% Volcanoplot
#!!! INDICATE ESSENTIAL GENES USING DIFFERENT MARKER
plt.figure(figsize=(19.0,9.0))#(27.0,3))
grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0,0])

# ax.scatter(x=volcano_df.fold_change, y=volcano_df.p_value, alpha=0.4, marker='.', c='k')
ax.scatter(x=volcano_df.loc[volcano_df['significance'] == False, 'fold_change'], y=volcano_df.loc[volcano_df['significance'] == False, 'p_value'], alpha=0.4, marker='.', c='k', label='p-value > {}'.format(significance_threshold))
ax.scatter(x=volcano_df.loc[volcano_df['significance'] == True, 'fold_change'], y=volcano_df.loc[volcano_df['significance'] == True, 'p_value'], alpha=0.4, marker='.', c='r', label='p-value < {}'.format(significance_threshold))
ax.grid(True, which='major', axis='both', alpha=0.4)
ax.set_xlabel('Log2 FC')
ax.set_ylabel('-1*Log10 p-value')
ax.set_title(variable)
ax.legend()
if not track_gene == "":
    ax.annotate(volcano_df.iloc[track_gene_index,:]['gene_names'], (volcano_df.iloc[track_gene_index,:]['fold_change'], volcano_df.iloc[track_gene_index,:]['p_value']),
                size=16, c='b')
    ax.scatter(x=volcano_df.iloc[track_gene_index,:]['fold_change'] ,y=volcano_df.iloc[track_gene_index,:]['p_value'], marker='X', c='b')

del (ax, grid)



#%%THIS IS AN ALTERNATIVE APPROACH FOR THE ABOVE CALCULATION.
### TEST INDEPENDENT T-TEST
### https://www.statisticshowto.com/independent-samples-t-test/

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

# print("t-value according to calculation:", t)
# print("t-value according to scipy:", stats.ttest_ind(test1,test2))












