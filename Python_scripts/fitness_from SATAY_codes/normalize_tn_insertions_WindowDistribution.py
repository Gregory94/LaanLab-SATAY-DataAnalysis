# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 10:51:28 2020

@author: gregoryvanbeek
"""

#%%
import os, sys
import numpy as np
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from genomicfeatures_dataframe_with_normalization import dna_features
from chromosome_and_gene_positions import chromosome_position


#%% LOAD FILES
gff_file = os.path.join(file_dirname,'..','Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
wig_file_path = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
pergene_insertions_file_path = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"


#%% USER SETTINGS
N_regions_perchrom_dict = {}
N_windows = 20000
chrom_list = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']

for chrom in chrom_list:
    print(chrom)
    chrom = chrom.upper()
    len_chr = chromosome_position(gff_file)[0].get(chrom)
    
    #%%
    dna_df2 = dna_features(region = chrom,
                 wig_file = wig_file_path,
                 pergene_insertions_file = pergene_insertions_file_path,
                 plotting=False, verbose=False)
    
    
    #%% DETERMINE START AND END POSITION OF THE REGIONS (WINDOWS)
    N_regions_dict = {}
    N_bp_dict = {}
    for N in [round(len_chr/N_windows)]:#[1,2,5,10,20,50,100]:
    
        regions_start = np.linspace(0, len_chr, N+1, dtype=int).tolist() #400 is approx 1000 bp windows. For a single region, use number of regions=2.
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
        
        N_bp_noncodingdna_list = [[]]
        N_bp_codingdna_list = [[]]
        N_bp_essgene_list = [[]]
        N_bp_nonessgene_list = [[]]
        for i in range(len(regions_list)):
            for f in dna_df2.itertuples(): #regions_list[i][0] is start position current region, region_list[i][1] is end position current region
                if regions_list[i][0] < f[6][0] < regions_list[i][1]: #f[4] contains and start and end position
                    if f[4] == None: #f[2] contains feature type
                        N_reads_noncodingdna_list[i].append(f[11]) #f[9] contains number of reads per bp
                        N_bp_noncodingdna_list[i].append(f[7]) #f[5] contains number of bp for this region
                    elif f[4].startswith("Gene") and f[5] == True: #f[3] contains whether a feature is essential or not
                        N_reads_essgene_list[i].append(f[11])
                        N_bp_essgene_list[i].append(f[7])
                    elif f[4].startswith("Gene") and f[5] == False:
                        N_reads_nonessgene_list[i].append(f[11])
                        N_bp_nonessgene_list[i].append(f[7])
                    else:
                        N_reads_codingdna_list[i].append(f[11])
                        N_bp_codingdna_list[i].append(f[7])
                        
                elif f[6][0] > regions_list[i][1] and i != len(regions_list):
                    N_reads_noncodingdna_list.append([])
                    N_reads_essgene_list.append([])
                    N_reads_nonessgene_list.append([])
                    N_reads_codingdna_list.append([])
                    
                    N_bp_noncodingdna_list.append([])
                    N_bp_essgene_list.append([])
                    N_bp_nonessgene_list.append([])
                    N_bp_codingdna_list.append([])
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
        bp = []
        for i in range(len(regions_list)):
            for j in N_reads_noncodingdna_list[i]:
                values.append(j)
                names.append('Noncoding_dna')
                regions.append(regions_list[i][0])
            for j in N_bp_noncodingdna_list[i]:
                bp.append(j)
               
            for j in N_reads_codingdna_list[i]:
                values.append(j)
                names.append('Coding_dna')
                regions.append(regions_list[i][0])
            for j in N_bp_codingdna_list[i]:
                bp.append(j)
        
            for j in N_reads_essgene_list[i]:
                values.append(j)
                names.append('Essential_gene')
                regions.append(regions_list[i][0])
            for j in N_bp_essgene_list[i]:
                bp.append(j)
        
            for j in N_reads_nonessgene_list[i]:
                values.append(j)
                names.append('Nonessential_gene')
                regions.append(regions_list[i][0])
            for j in N_bp_nonessgene_list[i]:
                bp.append(j)
        
        readsperbp_dict = {"Start_bp_per_region": regions,
                           "Genomic_type": names,
                           "N_basepairs": bp,
                           "Reads_per_bp": values}
        readsperbp_df = pd.DataFrame(readsperbp_dict, columns = [column_name for column_name in readsperbp_dict])
    #del (N_reads_noncodingdna_list, N_reads_codingdna_list, N_reads_essgene_list, N_reads_nonessgene_list, readsperbp_dict)
    
    
    
    #%% PRINT STATISTICS
        N_regions_list = []
        N_bp_list = []
    
        noncodingdna_bp = 0 #NUMBER OF NONCODING BASEPAIRS WITHIN THE WINDOWS
        readsperbp_noncoding_list = []
        noncodingdna_N = 0 #NUMBER OF NONCODING REGIONS WITHIN THE WINDOWS
        codingdna_bp = 0
        essgene_bp = 0
        nonessgene_bp = 0
        txt_box = 'Number of basepairs in each box\n'
        for reg in regions_list:
            for i in readsperbp_df.itertuples():
                if i[1] == reg[0]:
                    if i[2] == "Noncoding_dna":
                        noncodingdna_bp += i[3]
                        readsperbp_noncoding_list.append(i[4])
                        noncodingdna_N += 1
                    elif i[2] == "Coding_dna":
                        codingdna_bp += i[3]
                    elif i[2] == "Essential_gene":
                        essgene_bp += i[3]
                    else:
                        nonessgene_bp += i[3]
                elif i[1] > reg[0]:
    #                print('Number of noncoding regions in window %i = %i. Std in number of reads = %.2f.' % (reg[0], noncodingdna_N, np.std(readsperbp_noncoding_list)))
                    readsperbp_noncoding_list = []
                    txt_box = txt_box + 'N_' + str(reg[0]) + ' = ' + str(noncodingdna_bp) + ', '
                    txt_box = txt_box + str(codingdna_bp) + ', '
                    txt_box = txt_box + str(essgene_bp) + ', '
                    txt_box = txt_box + str(nonessgene_bp) + '\n'
                    
                    N_regions_list.append(noncodingdna_N)
                    N_bp_list.append(noncodingdna_bp)
                    
                    noncodingdna_bp = 0
                    noncodingdna_N = 0
                    codingdna_bp = 0
                    essgene_bp = 0
                    nonessgene_bp = 0
                    break
    #    print('Number of noncoding regions in window %i = %i. Std in number of reads = %.2f.' % (reg[0], noncodingdna_N, np.std(readsperbp_noncoding_list)))
        txt_box = txt_box + 'N_' + str(reg[0]) + ' = ' + str(noncodingdna_bp) + ', '
        txt_box = txt_box + str(codingdna_bp) + ', '
        txt_box = txt_box + str(essgene_bp) + ', '
        txt_box = txt_box + str(nonessgene_bp) + '\n'
        N_regions_list.append(noncodingdna_N)
        N_bp_list.append(noncodingdna_bp)
    
    
        N_regions_dict[N] = N_regions_list
        N_bp_dict[N] = N_bp_list
    
    
#    print(N_regions_dict)
#    print(N_bp_dict)
    
    N_regions_perchrom_dict[chrom] = N_regions_list
N_regions_perchrom_dict[chrom] = N_regions_list
print(N_regions_perchrom_dict)
#%% TEST PLOT
#plt.figure(figsize=(19,9))
#grid = plt.GridSpec(1, len(regions_list), wspace=0.0, hspace=0.0)
#
##df1 = readsperbp_df.stack().apply(pd.Series).stack().astype(float).rename_axis(['Region', 'dna', None]).rename('Reads_per_bp').reset_index(['Region','dna']).reset_index(drop=True)
##df1 = df1[df1.dna != 'Region']
##df1.head() #https://stackoverflow.com/questions/39344167/grouped-boxplot-with-seaborn
#
#ax = sns.boxplot(x='Start_bp_per_region', y='Reads_per_bp', hue='Genomic_type', data=readsperbp_df, hue_order=["Noncoding_dna", "Coding_dna", "Nonessential_gene", "Essential_gene"], showfliers=True, notch=False)
#ax.grid(axis='y', which='both')
#plt.xticks(rotation=90)
#plt.title('Chromosome ' + chrom)
#
#
##if N > 15:
##    ax.text(0.01, 0.5, txt_box, transform=ax.transAxes, fontsize=10, bbox=dict(boxstyle='round', facecolor='#c6e4e3', alpha=0.9))
##elif N < 5:
##    ax.text(0.01, 0.8, txt_box, transform=ax.transAxes, fontsize=10, bbox=dict(boxstyle='round', facecolor='#c6e4e3', alpha=0.9))
##else:
##    ax.text(0.01, 0.7, txt_box, transform=ax.transAxes, fontsize=10, bbox=dict(boxstyle='round', facecolor='#c6e4e3', alpha=0.9))
#                                                                                   
#
#
#plt.show()

#%% PLOT BOXPLOT FOR ALL CHROMOSOMES USING THE SAME WINDOW SIZE.

labels, data = [*zip(*N_regions_perchrom_dict.items())]

textsize = 20

plt.figure(figsize=(19,9))
grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.01)


ax = plt.subplot(grid[0,0])

c1 = (0.0,0.0,0.0,1.0)
c2 = '#F20064'
bp = ax.boxplot(data, notch=False, sym='.', patch_artist=True,
                boxprops=dict(facecolor='white', color=c1, linewidth=1),
                medianprops=dict(linewidth=2, color=c2),
                flierprops=dict(color=c2, markeredgecolor=c2, markersize=10),
                capprops=dict(color=c1, linewidth=1),
                whiskerprops=dict(color=c1, linewidth=1))

plt.xticks(list(range(1, len(labels) + 1)), chrom_list)
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

ax.set_xlabel('Chromosome number', fontsize=textsize)
ax.set_ylabel('Numer of neutral regions in each window', fontsize=textsize)
ax.set_title('Window length approximately %i' % N_windows, fontsize=textsize)

ax.tick_params(labelsize=textsize)

plt.show()











