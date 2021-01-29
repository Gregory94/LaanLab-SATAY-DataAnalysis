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
    file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt2_testfolder\wt2_matlab_tnmapping_20210127\E-MTAB-4885.WT2.bam.bed"
#    file=r"X:\tnw\BN\LL\Shared\Gregory\datasets\wt1_testfolder\E-MTAB-4885.WT1.sorted.bam.bed"
#    file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\wt1_dataset_enzo\wt1_enzo_dataset_demultiplexed_singleend_sample2_trim1\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam.bed"
    read_list = library_stats(bedfile=file)