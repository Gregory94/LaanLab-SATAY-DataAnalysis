# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a test script file for translating the Matlab code from the Kornmann lab to Python.

To save  numpy array in text format:
    np.savetxt(os.path.join(dirname,'tncoordinatescopy_python.txt'),tncoordinatescopy_array,delimiter=',',fmt='%i')
"""

import os, sys
import numpy as np
import pysam
import timeit

dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from chromosome_and_gene_positions import chromosomename_roman_to_arabic, gene_position
from gene_names import gene_aliases
from essential_genes_names import list_known_essentials

#%% START TIMER
timer_start_complete = timeit.default_timer()

#%% LOADING FILES
path = os.path.join('/home', 'gregoryvanbeek', 'Documents', 'data_processing')
# filename = os.path.join('E-MTAB-4885.WT2.bam')
filename = os.path.join('E-MTAB-4885.WT1.bam')

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
ref_name_list = []
for i in range(17):
    ref_name = bamfile.get_reference_name(i)
    ref_tid_dict[ref_name] = bamfile.get_tid(ref_name)
    ref_name_list.append(ref_name)

#%% GET SEQUENCE LENGTHS OF ALL CHROMOSOMES
chr_length_dict = {} # 'I' | 230218
chr_summedlength_dict = {}
ll = 0
for key in ref_tid_dict:
    ref_length = bamfile.get_reference_length(key)
    chr_length_dict[key] = ref_length
    chr_summedlength_dict[key] = ll
    ll += ref_length
del (ref_length, ll)

#%% GET NUMBER OF MAPPED, UNMAPPED AND TOTAL AMOUNT OF READS PER CHROMOSOME
chr_mappedreads_dict = {} # 'I' | [mapped, unmapped, total reads]
total_reads = 0
for i in range(17):
    stats = bamfile.get_index_statistics()[i]
    chr_mappedreads_dict[stats[0]] = [stats[1], stats[2], stats[3]]
    total_reads += stats[3]


#%% GET ALL READS WITHIN A SPECIFIED GENOMIC REGION
# readnumb_array = [] #np.array([], dtype=int)
tnnumber_dict = {}

ll = 0 #Number of unique insertions in entire genome
# temp = ['I','II']
for kk in ref_name_list: # 'kk' is chromosome number in roman numerals
    read_counter = 0
    timer_start = timeit.default_timer()
    
    N_reads_kk = chr_mappedreads_dict[kk][2]
    start_array = np.empty(shape=(N_reads_kk), dtype=int)
    flag_array = np.empty(shape=(N_reads_kk), dtype=int)
    readlength_array = np.empty(shape=(N_reads_kk), dtype=int)

    print('Getting reads for chromosme ', kk ,' ...')
    for reads in bamfile.fetch(kk, 0, chr_length_dict[kk]):
        read = str(reads).split('\t')
        
        start_array[read_counter] = int(read[3]) + 1
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
    
    
    start2_sortindices = start2_array.argsort(kind='mergesort') #use mergesort for stable sorting
    start2_array = start2_array[start2_sortindices]
    flag2_array = flag2_array[start2_sortindices]
    
    del start2_sortindices
    
    
##% CREATE ARRAY OF START POSITION AND FLAGS OF ALL READS IN GENOME
    ref_tid_kk = int(ref_tid_dict[kk]+1)
    if ll == 0:
        tncoordinates_array = np.array([])
    
    mm = 0 # Number of unique reads per insertion
    jj = 1 # Number of unique reads in current chromosome (Number of transposons in current chromosome)
    temp_counter = 0
    for ii in range(1,len(start2_array)):
        if abs(start2_array[ii]-start2_array[ii-1]) <= 2 and flag2_array[ii] == flag2_array[ii-1]:
            mm += 1
            temp_counter += 1
        else:
            avg_start_pos = abs(round(np.mean(start2_array[ii-mm-1 : ii])))
            if tncoordinates_array.size == 0:
                tncoordinates_array = np.array([ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])])
                readnumb_list = [mm+1]
            else:
                tncoordinates_array = np.vstack((tncoordinates_array, [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])]))    
                readnumb_list.append(mm+1)
            mm = 0
            jj += 1
            ll += 1
        
        if ii == len(start2_array) - 1: #include last read
            avg_start_pos = abs(round(np.mean(start2_array[ii-mm-1 : ii])))
            tncoordinates_array = np.vstack((tncoordinates_array, [ref_tid_kk, int(avg_start_pos), int(flag2_array[ii-1])]))
            readnumb_list.append(mm+1)

    tnnumber_dict[kk] = jj
    
    del jj, start_array, flag_array, readlength_array, flag2_array, start2_array

    timer_stop = timeit.default_timer()
    print('Loop over reads chromosome ', kk, ' complete. Time = ', timer_stop-timer_start, 'seconds')

readnumb_array = np.array(readnumb_list)
del readnumb_list

tncoordinatescopy_array = np.array(tncoordinates_array, copy=True)
# np.savetxt(os.path.join(dirname,'tncoordinates_python.txt'),tncoordinates_array,delimiter=',',fmt='%i') #!!!
# np.savetxt(os.path.join(dirname,'readnumb_python.txt'),readnumb_array,delimiter=',',fmt='%i') #!!!

#%% GET LIST OF ALL GENES AND ALL ESSENTIAL GENES

files_path = os.path.join(dirname,'..','data_files')
chromosomenames_gff_dict = chromosomename_roman_to_arabic()[0] #These names should correspond to the chromosome names in the gff files

# GET POSITION GENES
gff_path = os.path.join(files_path,'Saccharomyces_cerevisiae.R64-1-1.99.gff3')
genecoordinates_dict = gene_position(gff_path) #contains all genes, essential and nonessential
# gff_gene_name, gff_gene_chr, gff_gene_start, gff_gene_end = gene_position(gff_path, get_dict=False)[0:4]

# GET ALL ANNOTATED ESSENTIAL GENES
essential_path1 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_1.txt')
essential_path2 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_2.txt')
essentialnames_list = list_known_essentials([essential_path1, essential_path2])

# GET ALIASES OF ALL GENES
names_path = os.path.join(files_path,'Yeast_Protein_Names.txt')
aliases_designation_dict = gene_aliases(names_path)[0]

# FOR ALL GENE IN GENECOORDINATES_DICT, CHECK IF THEY ARE ANNOTATED AS ESSENTIAL
essentialcoordinates_dict = {}
gene_counter = 0
for gene in genecoordinates_dict:
# for gene in gff_gene_name:
    if gene in essentialnames_list:
        essentialcoordinates_dict[gene] = genecoordinates_dict.get(gene).copy()
        # essentialcoordinates_dict[gene] = [gff_gene_chr[gene_counter], gff_gene_start[gene_counter], gff_gene_end[gene_counter]]
    else:
        gene_aliases_list = []
        for key, val in aliases_designation_dict.items():
            if gene == key or gene in val: #if gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
                gene_aliases_list.append(key)
                for aliases in aliases_designation_dict.get(key):
                    gene_aliases_list.append(aliases)

        for gene_alias in gene_aliases_list:
            if gene_alias in essentialnames_list:
                essentialcoordinates_dict[gene_alias] = genecoordinates_dict.get(gene).copy()
                # essentialcoordinates_dict[gene_alias] = [gff_gene_chr[gene_counter], gff_gene_start[gene_counter], gff_gene_end[gene_counter]]
                break
    gene_counter += 1

#%% CONCATENATE ALL CHROMOSOMES
ll = 0
for ii in range(1,len(ref_name_list)):
    ll += chr_length_dict[ref_name_list[ii-1]]
    aa = np.where(tncoordinatescopy_array[:,0] == ii + 1)
    tncoordinatescopy_array[aa,1] = tncoordinatescopy_array[aa,1] + ll

for key in genecoordinates_dict:
    genecoordinates_dict[key][1] = genecoordinates_dict.get(key)[1] + chr_summedlength_dict.get(genecoordinates_dict.get(key)[0])
    genecoordinates_dict[key][2] = genecoordinates_dict.get(key)[2] + chr_summedlength_dict.get(genecoordinates_dict.get(key)[0])

for key in essentialcoordinates_dict:
    essentialcoordinates_dict[key][1] = essentialcoordinates_dict.get(key)[1] + chr_summedlength_dict.get(essentialcoordinates_dict.get(key)[0])
    essentialcoordinates_dict[key][2] = essentialcoordinates_dict.get(key)[2] + chr_summedlength_dict.get(essentialcoordinates_dict.get(key)[0])

# test = []
# for key in genecoordinates_dict:
#     test.append(genecoordinates_dict.get(key)[1:3])
#np.savetxt(os.path.join(dirname,'essentialcoordinates_python.txt'),test,delimiter=',',fmt='%i')

#%% GET NUMBER OF TRANSPOSONS AND READS PER GENE
tnpergene_dict = {}
readpergene_dict = {}
readpergenecrude_dict = {}
for gene in genecoordinates_dict:
    xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= genecoordinates_dict.get(gene)[1], tncoordinatescopy_array[:,1] <= genecoordinates_dict.get(gene)[2]))
    tnpergene_dict[gene] = np.size(xx)
    readpergene_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0)
    readpergenecrude_dict[gene] = sum(readnumb_array[xx])


tnperessential_dict = {}
readperessential_dict = {}
readperessentialcrude_dict = {}
for gene in essentialcoordinates_dict:
    xx = np.where(np.logical_and(tncoordinatescopy_array[:,1] >= essentialcoordinates_dict.get(gene)[1], tncoordinatescopy_array[:,1] <= essentialcoordinates_dict.get(gene)[2]))
    readperessential_dict[gene] = sum(readnumb_array[xx]) - max(readnumb_array[xx], default=0)
    readperessentialcrude_dict[gene] = sum(readnumb_array[xx])

#%% CREATE BED FILE
bedfile = file+'.bed'

with open(bedfile, 'w') as f:
    f.write('track name=' + filename + ' useScore=1\n')
    coordinates_counter = 0
    for tn in tncoordinates_array:
        refname = [key for key, val in ref_tid_dict.items() if val == tn[0] - 1][0]
        if refname == 'Mito':
            refname = 'M'
        f.write('chr' + refname + ' ' + str(tn[1]) + ' ' + str(tn[1] + 1) + ' . ' + str(100+readnumb_array[coordinates_counter]*20) + '\n')
        coordinates_counter += 1


#CHECK NUMBER OF INSERTIONS FOR EACH CHROMOSOME (NOTE, THIS IS SLOW)
# test_list = []
# for i in range(1,18):
#     test_counter = 0
#     for tn in tncoordinates_array:
#         if tn[0] == i:
#             test_counter += 1
#     test_list.append([i,test_counter])
# test_list

#%% CREATE TEXT FILE WITH TRANSPOSONS AND READS PER GENE
pergenefile = file+'_pergene.txt'

with open(pergenefile, 'w') as f:
    f.write('Gene name\tNumber of transposons per gene\tNumber of reads per gene\n')
    for gene in tnpergene_dict:
        tnpergene = tnpergene_dict[gene]
        readpergene = readpergene_dict[gene]
        f.write(gene + '\t' + str(tnpergene) + '\t' + str(readpergene) + '\n')

#%% CREATE WIG FILE
wigfile = file+'.wig'

readnumbwig_array = readnumb_array.copy()

unique_index_array = np.array([], dtype=int) #=cc
N_uniques_perchr_list = []
ll = 0
for kk in ref_name_list:
    index = np.where(tncoordinates_array[:,0] == int(ref_tid_dict[kk]+1)) #get indices for current chromosome.
    unique_index = np.unique(tncoordinates_array[index][:,1], return_index=True)[1] #get all insertion locations (in tncoordinates, all rows, column 1)

    unique_index_array = np.append(unique_index_array, (unique_index+ll), axis=0)
    
    if ll == 0:
        ll += np.count_nonzero(tncoordinates_array[:,0] == int(ref_tid_dict[kk]+1))# + 1
    else:
        ll += np.count_nonzero(tncoordinates_array[:,0] == int(ref_tid_dict[kk]+1))
    N_uniques_perchr_list.append(ll) #total amount unique indices found untill current chromosome

# np.savetxt(os.path.join(dirname,'uniqueindex_python.txt'),unique_index_array,delimiter=',',fmt='%i')

duplicate_list = [] #=dd
ll = 0
index_last_unique_previous_chromosome = 0
for ii in N_uniques_perchr_list:
    index_last_unique = np.where(unique_index_array <= ii)[0][-1]
    for jj in range(ll,ii):
        if int(jj) not in unique_index_array[index_last_unique_previous_chromosome:index_last_unique]:
            if ll == 0:
                duplicate_list.append(jj)
            else:
                duplicate_list.append(jj)#-1)
    index_last_unique_previous_chromosome = index_last_unique
    ll = ii
# np.savetxt(os.path.join(dirname,'duplicateindex_python.txt'),duplicate_list,delimiter=',',fmt='%i')
# testposition_dict = {}
# testorientation_dict = {}
# for ii in duplicate_list:
#     testposition_dict[ii] = tncoordinates_array[ii][1] - tncoordinates_array[ii-1][1]
#     testorientation_dict[ii] = tncoordinates_array[ii][2] - tncoordinates_array[ii-1][2]


#SUM READNUMB VALUES AT INDEX IN DUPLICATE_LIST AND DUPLICATE_LIST-1    
for ii in duplicate_list:
    readnumbwig_array[ii-1] = readnumbwig_array[ii-1] + readnumbwig_array[ii]

readnumbwig_duplicatesremoved_array = np.delete(readnumbwig_array, duplicate_list, axis=0)
tncoordinateswig_duplicatesremoved_array = np.delete(tncoordinates_array, duplicate_list, axis=0)


# del duplicate_list, unique_index_array, N_uniques_perchr_list


with  open(wigfile, 'w') as f:
    f.write('track type=wiggle_0 ,maxheightPixels=60 name='+filename+'\n')
    for kk in ref_name_list:
        f.write('VariableStep chrom=chr' + kk + '\n')
        
        index = np.where(tncoordinateswig_duplicatesremoved_array[:,0] == int(ref_tid_dict[kk]+1)) #get indices for current chromosome.
        for ii in index[0]:
            f.write(str(tncoordinateswig_duplicatesremoved_array[ii][1]) + ' ' + str(readnumbwig_duplicatesremoved_array[ii]) + '\n')
#%% END TIMER
timer_stop_complete = timeit.default_timer()
print('Script took %.2f seconds to run' %(timer_stop_complete - timer_start_complete))

