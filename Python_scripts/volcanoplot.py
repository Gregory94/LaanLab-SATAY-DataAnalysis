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

Code for showing gene name when hovering over datapoint is based on:
    - https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib

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
# datapath_a = 
# filenames_a = []
# datapath_b = 
# filenames_b = []
datapath_a = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing"
filenames_a = ["WT-a_pergene.txt", "WT-b_pergene.txt"]
datapath_b = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing"
filenames_b = ["dnrp1-1-a_pergene.txt", "dnrp1-1-b_pergene.txt", "dnrp1-2-a_pergene.txt", "dnrp1-2-b_pergene.txt"]


variable = 'read_per_gene' #'read_per_gene' 'tn_per_gene', 'Nreadsperinsrt'
significance_threshold = 0.01 #set threshold above which p-values are regarded significant
normalize=True

trackgene_list = ['nrp1']
# trackgene_list = ['bsd2', 'mtf2', 'abd1', 'pmp3','yjr087w','rpa49','osw7','atg23','gef1','mec1','yor293c-a']
# trackgene_list = ['cpr1','def1', 'rrt1', 'ssc1', 'rsc3','cup2']
# trackgene_list = ['cdc42', 'bem1', 'bem3', 'bem2', 'nrp1', 'cdc24', 'cla4', 'ste20']

#%%
def volcano(path_a, filelist_a, path_b, filelist_b, variable='read_per_gene', significance_threshold=0.01, normalize=True, trackgene_list=[]):
    '''
    This creates a volcano plot that shows the fold change between two libraries and the corresponding p-values.
    Input:
        - path_a, path_b: paths to location of the datafiles for library a and library b
        - filelist_a, filelist_b: list of the names of the datafiles for library a and library b located in path_a and path_b respectively
        - variable: tn_per_gene, read_per_gene or Nreadsperinsrt (default='read_per_gene')
        - significance_threshold: Treshold value above which the fold change is regarded significant, only for plotting (default=0.01)
        - normalize: Whether to normalize variable. If set to True, each gene is normalized based on the total count in each dataset (i.e. each file in filelist_) (default=True)
        - trackgene_list: Enter a list of gene name(s) which will be highlighted in the plot (e.g. ['cdc42', 'nrp1']). If empty list, no gene will be tracked. (default=[])

    Output:
        - volcano_df: pandas dataframe containing:
            - gene_names
            - fold change
            - t statistic
            - p value
            - whether p value is above threshold
        - volcanoplot with the log2 fold change between the two libraries and the -log10 p-value.

    Fold change is determined by the mean of dataset b (experimental set) divided by the mean of dataset a (reference set).
    The datasets can be of different length.
    P-value is determined based on the student t-test (scipy.stats.ttest_ind).

    NOTE: the log-fold change is currently not determined for genes that have zero insertions in either the reference or experimental set.

    Dependencies:
        - numpy
        - scipy
        - matplotlib
        - python_modules/dataframe_from_pergene.py (https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/master/python_modules/dataframe_from_pergene.py)
    '''

    datafiles_list_a = []
    datafiles_list_b = []
    for files in filelist_a:
        datafile = os.path.join(path_a, files)
        assert os.path.isfile(datafile), 'File not found at: %s' % datafile
        datafiles_list_a.append(datafile)
    for files in filelist_b:
        datafile = os.path.join(path_b, files)
        assert os.path.isfile(datafile), 'File not found at: %s' % datafile
        datafiles_list_b.append(datafile)

    del (files, datafile, path_a, path_b, filelist_a, filelist_b)



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
    # if not trackgene_list == "":
    #     trackgene_list = trackgene_list.upper()
    #     print('')
    #     print("TRACKING VARIABLE %s" % trackgene_list)
    #     trackgene_list_index = tnread_gene_a.loc[tnread_gene_a['gene_names'] == trackgene_list].index[0]
    #     print("tnread_gene_a:", tnread_gene_a.loc[tnread_gene_a['gene_names'] == trackgene_list])
    #     print("tnread_gene_b:", tnread_gene_b.loc[tnread_gene_b['gene_names'] == trackgene_list])


    del (datafile_a, datafile_b, count, tnread_gene_b)
    if trackgene_list == []:
        del tnread_gene_a


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

        #Determine fold change per gene between libraries
        if np.mean(variable_b_array[count]) == 0 and np.mean(variable_a_array[count]) == 0:
            fc_list[count] = 0
        # elif np.mean(variable_b_array[count]) == 0 and np.mean(variable_a_array[count]) != 0:
        #     fc_list[count] = np.log2(0.01 / np.mean(variable_a_array[count]))
        # elif np.mean(variable_b_array[count]) != 0 and np.mean(variable_a_array[count]) == 0:
        #     fc_list[count] = np.log2(np.mean(variable_b_array[count]) / 0.01)
        elif np.mean(variable_b_array[count]) == 0:
            if variable == 'tn_per_gene':
                fc_list[count] = np.log2(5 / np.mean(variable_a_array[count]))
            elif variable == 'read_per_gene':
                fc_list[count] = np.log2(25 / np.mean(variable_a_array[count]))
        # elif np.mean(variable_b_array[count]) == 0 or np.mean(variable_a_array[count]) == 0:
            # fc_list[count] = np.nan
            # fc_list[count] = np.log2(max(np.mean(variable_a_array[count]), np.mean(variable_b_array[count])))
        elif np.mean(variable_a_array[count]) == 0:
            if variable == 'tn_per_gene':
                fc_list[count] = np.log2(np.mean(variable_b_array[count]) / 5)
            elif variable == 'read_per_gene':
                fc_list[count] = np.log2(np.mean(variable_b_array[count]) / 25)
        else:
            fc_list[count] = np.log2(np.mean(variable_b_array[count]) / np.mean(variable_a_array[count]))
            # fc_list[count] = np.mean(variable_b_array[count]) - np.mean(variable_a_array[count])

        #Take sum of number of insertions per library 
        # if sum(variable_a_array[count]) == 0 and sum(variable_b_array[count]) == 0:
        #     fc_list[count] = 0
        # elif sum(variable_a_array[count]) == 0 or sum(variable_b_array[count]) == 0:
        #     fc_list[count] = np.log2(max(sum(variable_a_array[count]) / norm_a, sum(variable_b_array[count]) / norm_b))
        # else:
        #     fc_list[count] = np.log2((sum(variable_a_array[count]) / norm_a) / (sum(variable_b_array[count]) / norm_b))


        ### printing specific genes
        # if not trackgene_list == "" and count == trackgene_list_index:
        #     print("variable_a_array:", variable_a_array[count])
        #     print("variable_b_index:", variable_b_array[count])
        #     print("t-test value:", ttest_val)
        #     print("mean a:", np.mean(variable_a_array[count]))
        #     print("mean b:", np.mean(variable_b_array[count]))
        #     print("fold change:", fc_list[count])


    volcano_df['fold_change'] = fc_list
    volcano_df['t_statistic'] = ttest_tval_list
    volcano_df['p_value'] = ttest_pval_list
    volcano_df['significance'] = signif_thres_list

    del(count, val, ttest_val, ttest_tval_list, ttest_pval_list, fc_list, signif_thres_list)
    if normalize == True:
        del (norm_a, norm_b)


#%% Volcanoplot
    fig = plt.figure(figsize=(19.0,9.0))#(27.0,3))
    grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
    ax = plt.subplot(grid[0,0])

    # ax.scatter(x=volcano_df.loc[volcano_df['significance'] == False, 'fold_change'], y=volcano_df.loc[volcano_df['significance'] == False, 'p_value'], alpha=0.4, marker='.', c='k', label='p-value > {}'.format(significance_threshold))
    # sc = ax.scatter(x=volcano_df.loc[volcano_df['significance'] == True, 'fold_change'], y=volcano_df.loc[volcano_df['significance'] == True, 'p_value'], alpha=0.4, marker='.', c='r', label='p-value < {}'.format(significance_threshold))
    colors = {False:'black', True:'red'}
    sc = ax.scatter(x=volcano_df['fold_change'], y=volcano_df['p_value'], alpha=0.4, marker='.', c=volcano_df['significance'].apply(lambda x:colors[x]))
    ax.grid(True, which='major', axis='both', alpha=0.4)
    ax.set_xlabel('Log2 FC')
    ax.set_ylabel('-1*Log10 p-value')
    ax.set_title(variable)
    ax.scatter(x=[],y=[],marker='.',color='black', label='p-value > {}'.format(significance_threshold)) #set empty scatterplot for legend
    ax.scatter(x=[],y=[],marker='.',color='red', label='p-value < {}'.format(significance_threshold)) #set empty scatterplot for legend
    ax.legend()
    if not trackgene_list == []:
        genenames_array = volcano_df['gene_names'].to_numpy()
        for trackgene in trackgene_list:
            trackgene = trackgene.upper()
            if trackgene in genenames_array:
                trackgene_index = tnread_gene_a.loc[tnread_gene_a['gene_names'] == trackgene].index[0]
                trackgene_annot = ax.annotate(volcano_df.iloc[trackgene_index,:]['gene_names'], (volcano_df.iloc[trackgene_index,:]['fold_change'], volcano_df.iloc[trackgene_index,:]['p_value']),
                            size=10, c='green', bbox=dict(boxstyle="round", fc="w"))
                trackgene_annot.get_bbox_patch().set_alpha(0.6)
            else:
                print('WARNING: %s not found' % trackgene)
        del (tnread_gene_a, genenames_array)


    names = volcano_df['gene_names'].to_numpy()
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        # text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), 
        #                         " ".join([names[n] for n in ind["ind"]]))
        text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
        annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        # annot.get_bbox_patch().set_alpha(0.4)


    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)


#%% return function
    return(volcano_df)


#%%
if __name__ == '__main__':
    volcano_df = volcano(path_a=datapath_a, filelist_a=filenames_a,
            path_b=datapath_b, filelist_b=filenames_b,
            variable=variable,
            significance_threshold=significance_threshold,
            normalize=normalize,
            trackgene_list=trackgene_list)



#%% TEST INDEPENDENT T-TEST
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












