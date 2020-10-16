# -*- coding: utf-8 -*-
"""
"""

import os

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
            chrom_line = 'variablestep'
            for line in lines[1:]:
                if not line.lower().startswith(chrom_line):
                    line_split = line.strip('\n').split(' ')
                    total_reads += int(line_split[1])

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
            return(total_reads)



        elif extension[1] == '.bed':

            with open(file,'r') as f:
                lines = f.readlines()

            total_reads = 0
            for line in lines[1:]:
                line_split = line.strip('\n').split(' ')
                total_reads += (int(line_split[4]) - 100)/20

            if verbose == True:
                print('Number of mapped reads found = %i' % total_reads)
            return(total_reads)


        else:
            print('WARNING: No valid file extension. Files must have either of the following extension; .wig or .bed')



#%%
if __name__ == '__main__':
    total_mapped_reads(file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig", verbose=True)
#    total_mapped_reads(file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed", verbose=True)