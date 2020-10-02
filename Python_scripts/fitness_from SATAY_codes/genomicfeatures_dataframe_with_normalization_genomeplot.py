# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:35:59 2020

@author: gregoryvanbeek
"""
import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))

from genomicfeatures_dataframe_with_normalization import dna_features

sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from chromosome_and_gene_positions import chromosome_position


#%% SET PATHS TO FILES
region = "ii"
wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
pergene_insertions_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"
variable="reads"
normalize=True
normalization_window_size = 10000
plotting=False
verbose=False


#%% OPEN DNA_DF2

dna_df2 = dna_features(region, wig_file, pergene_insertions_file, variable, normalize, normalization_window_size, plotting, verbose)


Nreadsperbp_list = dna_df2['Nreadsperbp'].tolist()
feature_position_list = dna_df2['position'].tolist()

#%% BINNING OF THE READS
chr_length_dict = chromosome_position()[0]
chr_length = chr_length_dict.get(region.upper())







