# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:57:39 2020

@author: gregoryvanbeek
"""

import matplotlib.pyplot as plt
from chromosome_and_gene_positions import chromosome_position

#%% CHECK DISTRIBUTION OF READ COUNT

wig_file = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\satay_analysis_testdata\Output_Processing_WT1_KornmannLab\ERR1533147_trimmed.sorted.bam.wig"

chrom_list = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
chr_length_dict = chromosome_position(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Saccharomyces_cerevisiae.R64-1-1.99.gff3")[0]

with open(wig_file,'r') as f:
    lines = f.readlines()


read_list = [0] * sum([chr_len for chr_name, chr_len in chr_length_dict.items()])
c = 0
chrom_sum_len = 0
for line in lines[2:]:
    if c < 16:
        if not line.startswith('Variable'):
            line_split = line.strip('\n').split(' ')
            read_list[int(line_split[0]) + chrom_sum_len] = int(line_split[1])
        else:
            chrom_sum_len += chr_length_dict.get(chrom_list[c])
            c += 1

plt.hist(read_list, bins=range(min(read_list), max(read_list) + 100, 100))
