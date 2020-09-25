# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 10:51:28 2020

@author: gregoryvanbeek
"""

#%%
import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from normalize_tn_insertions import dna_features
from chromosome_and_gene_positions import chromosome_position

#%%
chrom = 'V'

gff_file = os.path.join(file_dirname,'..','Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
wig_file_path = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
pergene_insertions_file_path = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"
#%%
dna_df2 = dna_features(region = chrom,
             wig_file = wig_file_path,
             pergene_insertions_file = pergene_insertions_file_path,
             plotting=False)

len_chr = chromosome_position(gff_file)[0].get(chrom)

#%% DETERMINE START AND END POSITION OF THE REGIONS (WINDOWS)
regions_start = np.linspace(0, len_chr, 10, dtype=int).tolist() #400 is approx 1000 bp windows.
regions_start.pop() #remove last element
regions_end = []
for r in regions_start[1:]:
    regions_end.append(r-1)
regions_end.append(len_chr)
regions_list = []
for r in range(len(regions_start)):
    regions_list.append([regions_start[r], regions_end[r]])


print("%i regions with length = %i" % (len(regions_list), (regions_list[0][1] - regions_list[0][0])))
del (r, regions_start, regions_end)


#%%
N_reads_noncodingdna_list = [[]]
N_reads_codingdna_list = [[]]
N_reads_essgene_list = [[]]
N_reads_nonessgene_list = [[]]
for i in range(len(regions_list)):
    for f in dna_df2.itertuples():
        if regions_list[i][0] < f[4][0] < regions_list[i][1]: #regions_list[i][0] is start position current region, region_list[i][1] is end position current region
            if f[2] == None:
                N_reads_noncodingdna_list[i].append(f[9])
            elif f[2].startswith("Gene") and f[3] == True:
                N_reads_essgene_list[i].append(f[9])
            elif f[2].startswith("Gene") and f[3] == False:
                N_reads_nonessgene_list[i].append(f[9])
            else:
                N_reads_codingdna_list[i].append(f[9])
        elif f[4][0] > regions_list[i][1] and i != len(regions_list):
            N_reads_noncodingdna_list.append([])
            N_reads_essgene_list.append([])
            N_reads_nonessgene_list.append([])
            N_reads_codingdna_list.append([])
            break

del (i, f)
#%% CREATE DATAFRAME

#readsperbp_dict = {"Region": [i[0] for i in regions_list],
#                   "Noncodingdna": N_reads_noncodingdna_list,
#                   "Essentialgenes": N_reads_essgene_list,
#                   "Nonessentialgenes": N_reads_nonessgene_list,
#                   "Codingdna": N_reads_codingdna_list}

#readsperbp_df = pd.DataFrame(readsperbp_dict, columns = [column_name for column_name in readsperbp_dict])


values = []
names = []
regions = []
for i in range(len(regions_list)):
    for j in N_reads_noncodingdna_list[i]:
        values.append(j)
        names.append('Noncoding_dna')
        regions.append(regions_list[i][0])
    for j in N_reads_codingdna_list[i]:
        values.append(j)
        names.append('Coding_dna')
        regions.append(regions_list[i][0])
    for j in N_reads_essgene_list[i]:
        values.append(j)
        names.append('Essential_gene')
        regions.append(regions_list[i][0])
    for j in N_reads_nonessgene_list[i]:
        values.append(j)
        names.append('Nonessential_gene')
        regions.append(regions_list[i][0])

readsperbp_dict = {"Region": regions,
                   "Genomic_type": names,
                   "Reads_per_bp": values}
readsperbp_df = pd.DataFrame(readsperbp_dict, columns = [column_name for column_name in readsperbp_dict])
#del (N_reads_noncodingdna_list, N_reads_codingdna_list, N_reads_essgene_list, N_reads_nonessgene_list, readsperbp_dict)


#%% TEST PLOT
plt.figure(figsize=(19,9))
grid = plt.GridSpec(1, len(regions_list), wspace=0.0, hspace=0.0)

#df1 = readsperbp_df.stack().apply(pd.Series).stack().astype(float).rename_axis(['Region', 'dna', None]).rename('Reads_per_bp').reset_index(['Region','dna']).reset_index(drop=True)
#df1 = df1[df1.dna != 'Region']
#df1.head() #https://stackoverflow.com/questions/39344167/grouped-boxplot-with-seaborn

sns.boxplot(x='Region', y='Reads_per_bp', hue='Genomic_type', data=readsperbp_df, showfliers=False)














