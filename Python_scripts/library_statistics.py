# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:16:31 2021

@author: gregoryvanbeek

This code inputs a .bed file and determines
- the number of reads mapped
- number of transposons mapped
- the median number of reads per transposon.

"""

import os
import statistics

def library_stats(bedfile=None):

    assert os.path.isfile(bedfile), 'Bed file not found at: %s' % bedfile #check if given bed file exists

    reads_list = []
    with open(bedfile,'r') as f:
        for line in f:
            l = line.strip('\n').split(' ')
            if not l == [''] and not line.startswith('track'):
                reads_list.append((int(l[-1])-100)/20)

    print('Numer of insertions = %i' % len(reads_list))
    print('Number of reads = %i' % sum(reads_list))
    print('Median number of reads = %i' % statistics.median(reads_list))


    return(reads_list)

#%%
if __name__ == '__main__':
#    file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\ERR1533147_trimmed.sorted.bam.bed"
    file=r"C:\Users\gregoryvanbeek\Desktop\test_Matlab_wt1\E-MTAB-4885.WT1.bam.bed"
    read_list = library_stats(bedfile=file)
    