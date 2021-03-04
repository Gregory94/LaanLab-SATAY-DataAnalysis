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

    avg_ins_list = []
    reads_list = []
    previous_ins = 0
    with open(bedfile,'r') as f:
        for line in f:
            l = line.strip('\n').split(' ')
            if not l == [''] and not line.startswith('track'):
                reads_list.append((int(l[-1])-100)/20)
                
                current_ins = int(l[1])
                if not current_ins < previous_ins:
                    avg_ins_list.append(current_ins - previous_ins)
                previous_ins = current_ins



    print(os.path.basename(bedfile))
    print('Numer of insertions = %i' % len(reads_list))
    print('Number of reads = %i' % sum(reads_list))
    print('Median number of reads = %i' % statistics.median(reads_list))
    print('Average insertion frequency is %i insertions per bp' % statistics.mean(avg_ins_list))


    return(reads_list, avg_ins_list)

#%%
if __name__ == '__main__':
#    file=r"\\?\X:\tnw\BN\LL\Shared\Gregory\datasets\dataset_enzo\wt1_enzo_dataset_demultiplexed_interleaved_sample1\wt1_enzo_dataset_demultiplexed_singleend_sample1_trim20210127\align_out\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_singleend_trimmed.sorted.bam.bed"
#    file=r"\\?\X:\tnw\BN\LL\Shared\Gregory\datasets\dataset_enzo\wt1_enzo_dataset_demultiplexed_interleaved_sample2\wt1_enzo_dataset_demultiplexed_singleend_sample2_trim20210122\align_out\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_singleend_trimmed.sorted.bam.bed"
#    file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed"
    file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\testing_site\wt1_testfolder\analysis_matlab_kornmanncode\ERR1533147_trimmed.sorted.bam.bed"
    read_list, avg_ins_list = library_stats(bedfile=file)