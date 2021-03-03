# -*- coding: utf-8 -*-
"""
"""

import os
import numpy as np

#%%
def total_mapped_reads(file, verbose=False):
    '''
    This function gets the total number of reads that are present in the dataset.
    Input is either a .bed-file or a .wig-file.
    Output is an integer representing the total mapped reads in the dataset.
    '''
    if not os.path.isfile(file):
        print('WARNING: Following path not found: ', file)
        return()

    else:

        extension = os.path.splitext(file)



        if extension[1] == '.wig':
            
            with open(file, 'r') as f:
                lines = f.readlines()

            total_reads = 0
            total_ins = 0
            reads_list = []
            chrom_line = 'variablestep'
            for line in lines[1:]:
                if not line.lower().startswith(chrom_line):
                    line_split = line.strip('\n').split(' ')
                    reads_list.append(int(line_split[1]))
                    total_reads += int(line_split[1])
                    total_ins += 1

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
                print('Number of insertions found = %i' % total_ins)
                print('Median read per transposon = %i' % np.median(reads_list))
            return(total_reads)



        elif extension[1] == '.bed':

            with open(file,'r') as f:
                lines = f.readlines()

            total_reads = 0
            total_ins = 0
            reads_list = []
            for line in lines[1:]:
                line_split = line.strip('\n').split(' ')
                reads_list.append((int(line_split[4]) - 100)/20)
                total_reads += (int(line_split[4]) - 100)/20
                total_ins += 1

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
                print('Number of insertions found = %i' % total_ins)
                print('Median read per transposon = %i' % np.median(reads_list))
            return(total_reads)


        else:
            print('WARNING: No valid file extension. Files must have either of the following extension; .wig or .bed')



#%%
if __name__ == '__main__':
    # total_mapped_reads(file=r"C:\Users\gregoryvanbeek\Desktop\test_matlab_benoit\test_matlab_benoit_wt1\E-MTAB-4885.WT1.bam.wig", verbose=True)
    # total_mapped_reads(file=r"C:\Users\gregoryvanbeek\Desktop\test_matlab_benoit\test_matlab_benoit_wt1\E-MTAB-4885.WT1.bam.bed", verbose=True)
    total_mapped_reads(file=r"C:\Users\gregoryvanbeek\Desktop\test_matlab_benoit\test_matlab_benoit_wt2\E-MTAB-4885.WT2.bam.bed", verbose=True)