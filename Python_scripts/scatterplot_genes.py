# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 11:14:44 2021

@author: gregoryvanbeek
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from genomicfeatures_dataframe import dna_features


dna_df = dna_features(region = 3,
             wig_file = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam.wig",
             pergene_insertions_file = r"V:\tnw\bn\ll\Shared\Gregory\Labmeetings\Labmeeting_presentations\Labmeeting20210209_dataset\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam_pergene_insertions.txt",
             variable="reads",
             plotting=False,
             savefigure=False,
             verbose=True)


read_all_df = dna_df[['Essentiality', 'Nreadsperinsrt_truncatedgene']]
read_gene_df = read_all_df[(read_all_df.Essentiality == True) | (read_all_df.Essentiality == False)]


x_lin = np.linspace(0, len(read_gene_df)-1, len(read_gene_df))


plt.figure(figsize=(19,9))
ax = sb.scatterplot(x=x_lin, y=read_gene_df.Nreadsperinsrt_truncatedgene, hue=read_gene_df.Essentiality)
