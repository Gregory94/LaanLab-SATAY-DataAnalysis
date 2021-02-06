# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:14:44 2021

@author: gregoryvanbeek
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from genomicfeatures_dataframe import dna_features


#wigfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam.wig"
#pergeneinsertionsfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"
wigfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam.wig"
pergeneinsertionsfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"

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


plt.figure(figsize=(19,9))
colorpalette = sns.diverging_palette(10, 170, s=90, l=50, n=2) #https://seaborn.pydata.org/generated/seaborn.diverging_palette.html#seaborn.diverging_palette
ax = sns.scatterplot(x=x_lin, y=read_gene_df.Nreadsperinsrt_truncatedgene, hue=read_gene_df.Essentiality, palette=colorpalette, alpha=0.5, marker='|')
ax.grid(linestyle='-', alpha=0.8)
ax.set_xlim(-1, max(x_lin)+1)
ax.set_ylim(-1, 100)
ax.set_xticklabels([])
ax.set_ylabel('Reads per insertion')
ax.set_xlabel('Genes')
