# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:35:59 2020

@author: gregoryvanbeek
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))

sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from chromosome_and_gene_positions import chromosome_position
from genomicfeatures_dataframe_with_normalization import dna_features

#%% 
def genome_normalization_plot(wig_file=None, pergene_insertions_file=None, variable='reads', normalize=True, normalization_window_size=20000, plotting=False, verbose=False):

    if wig_file == None:
        wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
    if pergene_insertions_file == None:
        pergene_insertions_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"
    
    
    #%% OPEN DNA_DF2
    feature_position_list = []
    Nreadsperbp_list = []
    Nreadsperbp_normalized_list = []
    
    
    chr_length_dict = chromosome_position()[0]
    chr_length = 0
    
    chrom_list = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
    
    for region in chrom_list:
#        print(region)
        dna_df2 = dna_features(region, wig_file, pergene_insertions_file, variable, normalize, normalization_window_size, plotting, verbose)
    
        feature_position_list = feature_position_list + [pos[0] + chr_length for pos in dna_df2['position'].tolist()]
        Nreadsperbp_list = Nreadsperbp_list + dna_df2['Nreadsperbp'].tolist()
        Nreadsperbp_normalized_list = Nreadsperbp_normalized_list + dna_df2['Nreadsperbp_normalized_byNCregions'].tolist()
    
        chr_length += chr_length_dict.get(region.upper())
    
    #%% BINNING OF THE READS
    N_bins = 10000
    bins = np.linspace(0, chr_length, round(chr_length/N_bins), dtype=int).tolist()
    
    b_start = bins[0]
    binned_reads_list = np.zeros(len(bins))
    binned_norm_reads_list = np.zeros(len(bins))
    i = 0
    for b_end in bins[1:]:
        for ind, pos in enumerate(feature_position_list):
            if b_start <= pos < b_end:
                binned_reads_list[i] = binned_reads_list[i] + Nreadsperbp_list[ind]
                binned_norm_reads_list[i] = binned_norm_reads_list[i] + Nreadsperbp_normalized_list[ind]
            elif pos > b_end:
                b_start = b_end
                break
        i += 1
    
    
    #%% SUMMED CHROMOSOME LENGTHS FOR PLOTTING GRID
    
    l_genome = 0
    chr_summedlength_dict = {}
    #!!! THE ORDER IS WRONG
    for chrom, length in chr_length_dict.items():
        chr_summedlength_dict[chrom] = l_genome
        l_genome += length
    
    
    
    
    
    
    
    #%% PLOTTING
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.01)
    
    
    ax1 = plt.subplot(grid[0,0])
    ax1.bar(bins, binned_reads_list, width=N_bins, color="#333333")
    ax1.grid(False)
    ax1.set_xlim(0,l_genome)
    ax1.set_ylim(0,250)
    ax1.set_ylabel("Reads not normalized [Absolute counts]")
    
    ax2 = plt.subplot(grid[1,0])
    ax2.bar(bins, binned_norm_reads_list, width=N_bins, color="#333333")
    ax2.grid(False)
    ax2.set_xlim(0,l_genome)
    ax2.set_ylim(0,250)
    ax2.set_ylabel("Reads normalized [A.U.]")
    
    for chrom in chr_summedlength_dict:
        ax1.axvline(x = chr_summedlength_dict.get(chrom), linestyle='-', color=(0.9,0.9,0.9,1.0))
        ax2.axvline(x = chr_summedlength_dict.get(chrom), linestyle='-', color=(0.9,0.9,0.9,1.0))
    
    
    
    
    middle_chr_position = []
    c1 = chr_summedlength_dict.get('I')
    for c in chr_summedlength_dict:
        if not c == 'I':
            c2 = chr_summedlength_dict.get(c)
            middle_chr_position.append(c1 + (c2 - c1)/2)
            c1 = c2
    c2 = l_genome
    middle_chr_position.append(c1 + (c2 - c1)/2)
    
    ax2.set_xticks(middle_chr_position)
    ax2.set_xticklabels(chrom_list)
    ax2.tick_params(axis='x', which='major', pad=10)


#%%
if __name__ == '__main__':
    genome_normalization_plot()
