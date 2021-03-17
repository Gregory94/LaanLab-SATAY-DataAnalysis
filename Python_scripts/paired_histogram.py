# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 19:26:59 2021

@author: gregoryvanbeek
"""

import os, sys
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


#%% INPUT
datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing\WT-a_pergene.txt"




#%% CREATE HISTOGRAM COMPARISON
def hist_compare(datafile):

    assert os.path.isfile(datafile), 'File not found at: %s' % datafile

    with open(datafile) as f:
        lines = f.readlines()[1:] #skip header

    genenames_list = [None]*len(lines)
    tnpergene_list = [None]*len(lines)
    readpergene_list = [None]*len(lines) 

    line_counter = 0
    for line in lines:
#        l = line.strip('\n').split(' ')
        l = re.split(' |\t', line.strip('\n'))

        genenames_list[line_counter] = l[0]
        tnpergene_list[line_counter] = int(l[1])
        readpergene_list[line_counter] = int(l[2])

        line_counter += 1

    del (line, l, line_counter, datafile)


#%% determine number of reads per insertion per gene
    readperinspergene_list = [None]*len(lines)
    for i in range(len(tnpergene_list)):
        if not tnpergene_list[i] == 0:
            readperinspergene_list[i] = readpergene_list[i] / tnpergene_list[i]
        else:
            readperinspergene_list[i] = 0

    del (i)


#%% determine essential genes
    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
    from essential_genes_names import list_known_essentials #import essential_genes_names from python modules directory

    known_essential_gene_list = list_known_essentials(input_files=[os.path.join(file_dirname,'..','Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                                                                   os.path.join(file_dirname,'..','Data_Files','Cerevisiae_EssentialGenes_List_2.txt')])

    geneessentiality_list = [None]*len(lines)
    for i in range(len(genenames_list)):
        if genenames_list[i] in known_essential_gene_list:
            geneessentiality_list[i] = True
        else:
            geneessentiality_list[i] = False

    del (lines, file_dirname, known_essential_gene_list, i)

#%% create dataframe
    read_gene_dict = {"gene_names": genenames_list,
                      "gene_essentiality": geneessentiality_list,
                      "tn_per_gene": tnpergene_list,
                      "read_per_gene": readpergene_list,
                      "Nreadsperinsrt": readperinspergene_list}

    read_gene_df = pd.DataFrame(read_gene_dict, columns = [column_name for column_name in read_gene_dict])

    del (read_gene_dict, genenames_list, geneessentiality_list, tnpergene_list, readpergene_list, readperinspergene_list)



#%% create plot
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)

    ax = plt.subplot(grid[0,0])
    sns.histplot(x=read_gene_df.loc[read_gene_df["gene_essentiality"] == False, "tn_per_gene"], color='r', binwidth=5)
    max_x = ax.get_xlim()
    ax.set_xlim(left=0)
    ax.grid(True)

    ax2 = plt.subplot(grid[1,0])    
    sns.histplot(x=read_gene_df.loc[read_gene_df["gene_essentiality"] == True, "tn_per_gene"], color='g', binwidth=5)
    ax2.invert_yaxis()
    ax2.set_xlim(0,max_x[1])
    ax2.grid(True)

    plt.show()

#%% return
    return (read_gene_df)



#%%
if __name__ == '__main__':
    read_gene_df = hist_compare(datafile)

