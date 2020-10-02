# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:43:12 2020

@author: gregoryvanbeek
"""
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
#%%
file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"

with open(file, 'r') as f:
    lines = f.readlines()

gene_len_list = []
Ninsrt_list = []
for line in lines[1:]:
    line_split = line.strip('\n').split('\t')

    gene_len_list.append(int(line_split[3]) - int(line_split[2]))
    Ninsrt_list.append(len(line_split[4].replace('[', '').replace(']','').split(', ')))






#%%
plt.figure(figsize=(19,9))
grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)
ax = plt.subplot(grid[0:19,0])

ax.scatter(x=gene_len_list, y=Ninsrt_list, alpha=0.5, c='k')