# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:14:44 2021

@author: gregoryvanbeek


This script contains two functions that do the same thing, but have a different approach:
    - essentialgene_plot()
    - essentialgene_plot_frompergenefile()

The essentialgene_plot() function uses the genomicfeatures_dataframe.py to acquire the information regarding gene essentiality and number of reads per insertion per gene.
The essentialgene_plot_frompergenefile() uses the .bam_pergene.txt file that is directly outputted from transposonmapping_satay.py.

The main difference between the two files is that in the .bam_pergene.txt file the insertion with the highest number of reads are discarded to reduce the noise.
(For this, see the `readpergene_dict` variable in transposonmapping_satay.py and this discussion https://groups.google.com/g/satayusers/c/uaTpKsmgU6Q/m/IIHfSi0iBgAJ)
Therefore the graphs between the essentialgene_plot() and essentialgene_plot_frompergenefile() look slightly different, but the main distribution of the data should look similar.
Also the essentialgene_plot_frompergenefile() function typically has more genes with zero reads per insertion because of this feature.

When both files are available, it is adviced to start with the essentialgene_plot_frompergenefile().
To have the full dataset without the feature of removing the highest read count, use essentialgene_plot().
"""

import os, sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


#%%
def essentialgene_plot(wigfile, pergeneinsertionsfile):
    '''
    This code creates a scatterplot of the number of reads per insertion per gene combined with a histogram.
    The genes are sorted based on the number of reads per insertion and are color coded based on the annotated essentiality in wild type.

    Input:
        - path to .wig file
        - path to _pergene_insertions.txt file as created by transposonmapping_satay.py

    Requirements:
        - essential_genes_names.py located in python_modules directory (the python_modules directory is expected to be located in the same directory as this script).
        - Cerevisiae_EssentialGenes_List_1.txt and Cerevisiae_EssentialGenes_List_2.txt, located in the Data_Files directory (the Data_Files directory is expected to be located in the parent directory of this script).
    '''


    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
    from genomicfeatures_dataframe import dna_features

#%%
    assert os.path.isfile(wigfile), 'Wig file not found at: %s' % wigfile
    assert os.path.isfile(pergeneinsertionsfile), 'pergeneinsertionsfile file not found at: %s' % pergeneinsertionsfile

    loop_counter = 1
    for chrom in range(1,17):
        dna_df = dna_features(region = chrom,
                     wig_file = wigfile,
                     pergene_insertions_file = pergeneinsertionsfile,
                     variable="reads",
                     plotting=False,
                     savefigure=False,
                     verbose=True)


        read_all_df = dna_df[['Standard_name', 'Essentiality', 'Nreadsperinsrt_truncatedgene']]
        if loop_counter == 1:
            read_gene_df = read_all_df[(read_all_df.Essentiality == True) | (read_all_df.Essentiality == False)]
        else:
            read_gene_df_currentloop = read_all_df[(read_all_df.Essentiality == True) | (read_all_df.Essentiality == False)]
            read_gene_df = read_gene_df.append(read_gene_df_currentloop, ignore_index=True)


        loop_counter += 1


    read_gene_df = read_gene_df.sort_values(by=['Nreadsperinsrt_truncatedgene'])


    x_lin = np.linspace(0, len(read_gene_df)-1, len(read_gene_df))


    #%%
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(1, 20, wspace=0.0, hspace=0.0)


    ax1 = plt.subplot(grid[0,0:15])
    colorpalette = sns.diverging_palette(10, 170, s=90, l=50, n=2) #https://seaborn.pydata.org/generated/seaborn.diverging_palette.html#seaborn.diverging_palette
    sns.scatterplot(x=x_lin, y=read_gene_df.Nreadsperinsrt_truncatedgene, hue=read_gene_df.Essentiality, palette=colorpalette, alpha=0.5, marker='|', legend=True)
    ax1.grid(linestyle='-', alpha=0.8)
    ax1.set_xlim(-1, max(x_lin)+1)
    ax1.set_ylim(-1, 100)
    ax1.set_xticklabels([])
    ax1.set_ylabel('Reads per insertion')
    ax1.set_xlabel('Genes')

    ax2 = plt.subplot(grid[0,15:19])
    colorpalette = sns.diverging_palette(170, 10, s=90, l=50, n=2)
    sns.histplot(data=read_gene_df, y="Nreadsperinsrt_truncatedgene", hue="Essentiality", hue_order=[True, False], palette=colorpalette, alpha=0.5)
    ax2.get_legend().remove()
    ax2.set_xlim(0, 500)
    ax2.set_ylim(-1, 100)
    ax2.set_yticklabels([])
    ax2.set_ylabel('')
    ax2.grid(linestyle='-', alpha=0.8)


#%%
    return(read_gene_df)


#%%
#
#
#
#
#
#
#
#
#%%
def essentialgene_plot_frompergenefile(datafile):
    '''
    This code creates a scatterplot of the number of reads per insertion per gene combined with a histogram.
    The genes are sorted based on the number of reads per insertion and are color coded based on the annotated essentiality in wild type.
    
    Input:
        - path to _pergene.txt (each line containing a gene with corresponding number of insertions and reads seperated either by a space or tab)

    Requirements:
        - essential_genes_names.py located in python_modules directory (the python_modules directory is expected to be located in the same directory as this script).
        - Cerevisiae_EssentialGenes_List_1.txt and Cerevisiae_EssentialGenes_List_2.txt, located in the Data_Files directory (the Data_Files directory is expected to be located in the parent directory of this script).
    '''

    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
    from essential_genes_names import list_known_essentials #import essential_genes_names from python modules directory


#%%
#THIS COMMENT IS FOR THE FOLLOWING SECTIONS:
#    - READ FILE
#    - DETERMINE NUMBER OF READS PER INSERTION PER GENE
#    - DETERMINE ESSENTIAL GENES
#    - CREATE DATAFRAME
#
#THIS CAN BE REPLACED BY THE FOLLOWNG CODE:
#    file_dirname = os.path.dirname(os.path.abspath('__file__'))
#    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
#    from dataframe_from_pergene import dataframe_from_pergenefile
#    
#    read_gene_df_a = dataframe_from_pergenefile(datafile_a)


#%% read file
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


#%% sort values
    read_gene_df = read_gene_df.sort_values(by=["Nreadsperinsrt"])

    x_lin = np.linspace(0, len(read_gene_df)-1, len(read_gene_df))


#%% plotting
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(1, 20, wspace=0.0, hspace=0.0)


    ax1 = plt.subplot(grid[0,0:15])
    colorpalette = sns.diverging_palette(10, 170, s=90, l=50, n=2) #https://seaborn.pydata.org/generated/seaborn.diverging_palette.html#seaborn.diverging_palette
    sns.scatterplot(x=x_lin, y=read_gene_df.Nreadsperinsrt, hue=read_gene_df.gene_essentiality, palette=colorpalette, alpha=0.5, marker='|', legend=True)
    ax1.grid(linestyle='-', alpha=0.8)
    ax1.set_xlim(-1, max(x_lin)+1)
    ax1.set_ylim(-1, 100)
    ax1.set_xticklabels([])
    ax1.set_ylabel('Reads per insertion')
    ax1.set_xlabel('Genes')

    ax2 = plt.subplot(grid[0,15:19])
    colorpalette = sns.diverging_palette(170, 10, s=90, l=50, n=2)
    sns.histplot(data=read_gene_df, y="Nreadsperinsrt", hue="gene_essentiality", hue_order=[True, False], palette=colorpalette, alpha=0.5, binwidth=1)
    ax2.get_legend().remove()
#    ax2.set_xlim(0, 500)
    ax2.set_ylim(-1, 100)
    ax2.set_yticklabels([])
    ax2.set_ylabel('')
    ax2.grid(linestyle='-', alpha=0.8)


#%%
    return(read_gene_df)


#%%
#
#
#
#
#
#
#
#
#%% example inputs for data files for essentialgene_plot()

### dBEM1dBEM2dBEM3dNRP1 enzo
#wigfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam.wig"
#pergeneinsertionsfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"

### WILD TYPE enzo
#wigfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample2_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam.wig"
#pergeneinsertionsfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample2_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"

### WILD TYPE leila
# wigfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\leila_dataset_wt_processing\WT_merged-techrep-a_techrep-b_processing\WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_clean.wig"
# pergeneinsertionsfile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\leila_dataset_wt_processing\WT_merged-techrep-a_techrep-b_processing\WT_merged-techrep-a_techrep-b_trimmed.sorted.bam_pergene_insertions.txt"
#%% example inputs for data files for essentialgene_plot_frompergenefile()

### WILD TYPE leila
datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing\WT-a_pergene.txt"
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\dataset_leila_wt_agnesprocessing\WT-b_pergene.txt"

### dNRP1 leila
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing\dnrp1-1-a_pergene.txt"
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing\dnrp1-1-b_pergene.txt"
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing\dnrp1-2-a_pergene.txt"
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_dnpr1\dataset_leila_dnrp1_agnesprocessing\dnrp1-2-b_pergene.txt"

### dBEM1dBEM2dBEM3dNRP1 enzo
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam_pergene.txt"

#### WILD TYPE enzo
#datafile = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample2_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam_pergene.txt"

#%%
if __name__ == '__main__':
    # read_gene_df = essentialgene_plot(wigfile, pergeneinsertionsfile)

    read_gene_df = essentialgene_plot_frompergenefile(datafile)



