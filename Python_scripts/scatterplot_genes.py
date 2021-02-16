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

#%%
wigfile = r"\\?\N:\tnw\BN\LL\Shared\Gregory\datasets\dataset_enzo\wt1_enzo_dataset_demultiplexed_interleaved_sample1\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim20210127\align_out\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam.wig"
pergeneinsertionsfile = r"\\?\N:\tnw\BN\LL\Shared\Gregory\datasets\dataset_enzo\wt1_enzo_dataset_demultiplexed_interleaved_sample1\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim20210127\align_out\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"
#wigfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam.wig"
#pergeneinsertionsfile = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt"


#%%
def essentialgene_plot(wigfile, pergeneinsertionsfile):

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
#    ax1.legend(loc='upper left')
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
if __name__ == '__main__':
    essentialgene_plot(wigfile, pergeneinsertionsfile)




