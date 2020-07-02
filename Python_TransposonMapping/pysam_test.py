# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a test script file for translating the Matlab code from the Kornmann lab to Python.
"""

import os, sys
import numpy as np
import pysam

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,file_dirname)
# from chromosome_and_gene_positions import gene_position

#%% LOADING FILES
path = os.path.join('/home', 'gregoryvanbeek', 'Documents', 'data_processing')
filename = os.path.join('E-MTAB-4885.WT2.bam')

file = os.path.join(path,filename)
print('Running path: ', file)

if os.path.isfile(file):
    print('File exists.')
elif os.path.exists(file):
    print('File does not exist, but path does exists.')
else:
    print('Path does not exist.')


#%% READ BAM FILE
bamfile = pysam.AlignmentFile(file, 'rb')


#%% GET NAMES OF ALL CHROMOSOMES AS STORED IN THE BAM FILE
ref_tid_dict = {} # 'I' | 0
for i in range(17):
    ref_name = bamfile.get_reference_name(i)
    ref_tid_dict[ref_name] = bamfile.get_tid(ref_name)


#%% GET SEQUENCE LENGTHS OF ALL CHROMOSOMES
chr_length_dict = {} # 'I' | 230218
for key in ref_tid_dict:
    chr_length_dict[key] = bamfile.get_reference_length(key)


#%% GET NUMBER OF MAPPED, UNMAPPED AND TOTAL AMOUNT OF READS PER CHROMOSOME
chr_mappedreads_dict = {} # 'I' | [mapped, unmapped, total reads]
total_reads = 0
for i in range(17):
    stats = bamfile.get_index_statistics()[i]
    chr_mappedreads_dict[stats[0]] = [stats[1], stats[2], stats[3]]
    total_reads += stats[3]


#%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION
start_array = np.empty(shape=(total_reads), dtype=int)
flag_array = np.empty(shape=(total_reads), dtype=int)
readlength_array = np.empty(shape=(total_reads), dtype=int)

readnumb_array = np.array([], dtype=int)
tnnumber_dict = {}

read_counter = 0
ll = 0 #Number of unique insertions in entire genome
temp = 'I'
for kk in temp:#ref_tid_dict: # 'kk' is chromosome number in roman numerals
    print('Gettings reads for chromosme ', key)
    for reads in bamfile.fetch(kk, 0, chr_length_dict[kk]):
        read = str(reads).split('\t')
        
        start_array[read_counter] = int(read[3]) +1
        flag_array[read_counter] = int(read[1])
        readlength_array[read_counter] = int(len(read[9]))

        read_counter += 1


##% CORRECT STARTING POSITION FOR READS WITH REVERSED ORIENTATION
    flag0coor_array = np.where(flag_array==0) #coordinates reads 5' -> 3'
    flag16coor_array = np.where(flag_array==16) # coordinates reads 3' -> 5'

    startdirect_array = start_array[flag0coor_array]
    flagdirect_array = flag_array[flag0coor_array]

    startindirect_array = start_array[flag16coor_array] + readlength_array[flag16coor_array]
    flagindirect_array = flag_array[flag16coor_array]

    start2_array = np.concatenate((startdirect_array, startindirect_array), axis=0)
    flag2_array = np.concatenate((flagdirect_array, flagindirect_array), axis=0)
    
    del flag0coor_array, flag16coor_array, startdirect_array, flagdirect_array, startindirect_array, flagindirect_array
    
    
    start2_sortindices = start2_array.argsort()
    start2_array = start2_array[start2_sortindices]
    flag2_array = flag2_array[start2_sortindices]
    
    del start2_sortindices
    
    
##% CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENEOME
    tncoordinates_array = np.array([ref_tid_dict[kk]+1, start2_array[0], flag2_array[0]])
    
    mm = 0 # Number of unique reads per insertion
    jj = 1 # Number of unique reads in current chromosome (Number of transposons in current chromosome)
    for ii in range(1,len(start2_array)):
        if start2_array[ii]-start2_array[ii-1] <= 2 and flag2_array[ii] == flag2_array[ii-1]:
            mm += 1
        else:
            avg_start_pos = abs(int(np.mean(start2_array[ii-mm-1 : ii]))) #!!!CHECK IF LAST VALUE IS CONSIDERED CORRECTLY!!!
            tncoordinates_array = np.vstack((tncoordinates_array), [kk, avg_start_pos, flag2_array[ii]])
            mm = 0
            jj += 1
            ll += 1

        readnumb_array = np.append(readnumb_array, mm+1)

    tnnumber_dict[kk] = jj
    
    del jj, start_array, flag_array, readlength_array, start2_array, flag2_array


