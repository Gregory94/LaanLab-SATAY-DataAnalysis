# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:25:54 2020

@author: gregoryvanbeek
"""

import os, sys
import numpy as np
from tabulate import tabulate
sys.path.insert(1,r'C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\python_modules')
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic
from chromosome_names_in_files import chromosome_name_bedfile

#%%
def textfile_dataset_compare(bed_files=None,text_file=False):
    '''Create a text file that includes statistics for all chromosomes and the entire genome.
    The input 'bed_file' must be a list that include either one or two paths to a bed_file that needs to be analyzed
    The input text_file must be a string that indicates a path to the file where the data is written to.
    If text_file is empty, it will create a file with the same name as the first bed_file with the extension '_DataCompare.txt' at the same location.
    '''

#%% USED FILES
    gff_file = r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3"

#%% SET FILE NAME FOR CREATED TEXT FILE
    if text_file == True:
        text_file = os.path.splitext(bed_files[0])[0] + '_DataCompare.txt'
    
    if text_file != False:
        t = open(text_file, 'w+')
        if len(bed_files) == 1:
            t.write('Statistical values for datafile: ' + bed_files[0].split('\\')[-1] + '\n')
        elif len(bed_files) == 2:
            t.write('Statistical values for datafile: ' + bed_files[0].split('\\')[-1] + ' and ' + bed_files[1].split('\\')[-1] + '\n')
        t.close()
#%% GET CHROMOSOME START AND END POSTIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)

#%% GET ROMAN ARABIC NUMERALS
    roman_to_arabic_dict = chromosomename_roman_to_arabic()[1]
    chromosome_romannames_list = []
    for roman in roman_to_arabic_dict:
        chromosome_romannames_list.append(roman)

#%% ALLOCATE LISTS FOR PRINTING
#!!!CHANGE NAMING OF THESE VARIABLES TO SOMETHING MORE LOGIC
    line1 = []
    line2 = []
    line3 = []
    line4 = []
    line5 = []
    line6 = []
    line7 = []
    line8 = []
    line9 = []
    line10 = []
    line11 = []

#%% OPEN BED FILE
    bed_file_counter = 0
    for bed_file in bed_files:
        with open(bed_file) as f:
            lines = f.readlines()

#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
        chrom_names_dict, chrom_start_index_dict, chrom_end_index_dict = chromosome_name_bedfile(lines)

#%% DETERMINE STATISTICS FOR INDIVIDUAL CHROMOSOMES AND PUT THE RESULTS IN A DICT
#        if chromosome != None:#IF A SPECIFIC CHROMOSOME IS REQUISTED, TAKE ONLY THAT ONE
#            chromosome = chromosome.upper()
#            chrom_loop = {}
#            chrom_loop[chromosome] = chrom_names_dict.get(chromosome)
#        else: #IF NO SPECIFIC CHROMOSOME IS REQUISTED, TAKE ALL CHROMOSOMES
        chrom_loop = chrom_names_dict
    
        bp_between_tn_insertions_dict = {}
        reads_per_tn_dict = {}
        for chrom in chrom_loop:
            tn_insertion_position_list = []
            reads_per_tn_list = []
            for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
                line = line.strip('\n').split()
                tn_insertion_position_list.append(int(line[1]))
                reads_per_tn_list.append((int(line[4])-100)/20)
            bp_between_tn_insertions = [abs(y-x) for x, y in zip(tn_insertion_position_list[:-1], tn_insertion_position_list[1:])]
            bp_between_tn_insertions.insert(0,tn_insertion_position_list[0]) #ADD START OF GENE (bp=0)
            bp_between_tn_insertions.append(chr_length_dict.get(chrom) - tn_insertion_position_list[-1]) #ADD END OF GENE (bp=INDEX LAST TN - GENE LENGTH)
            bp_between_tn_insertions_dict[chrom] = bp_between_tn_insertions
            reads_per_tn_dict[chrom] = reads_per_tn_list
    
            tn_insertion_meanfrequency = np.nanmean(bp_between_tn_insertions)
            tn_insertion_25percentilefrequency = np.percentile(bp_between_tn_insertions,25)
            tn_insertion_medianfrequency = np.nanmedian(bp_between_tn_insertions)
            tn_insertion_75percentilefrequency = np.percentile(bp_between_tn_insertions,75)


            if bed_file_counter == 0:
                print('Print information chromosome ' + chrom + ' with length ' + str(chr_length_dict.get(chrom)))
                line1.append([chrom, 'Number of transposon insertions', len(reads_per_tn_list), ''])
                line2.append([chrom, 'Coverage percent', len(tn_insertion_position_list)/chr_length_dict.get(chrom)*100, ''])
                
                line3.append([chrom, 'Mean distance between insertions', tn_insertion_meanfrequency, ''])
                line4.append([chrom, 'Median distance between insertions', tn_insertion_medianfrequency, ''])
                line5.append([chrom, '25th percentile distance between insertions', tn_insertion_25percentilefrequency, ''])
                line6.append([chrom, '75th percentile distance between insertions', tn_insertion_75percentilefrequency, ''])
                
                line7.append([chrom, 'Largest area devoid of transposons', max(bp_between_tn_insertions), ''])
                
                line8.append([chrom, 'Mean number of reads per transposon', np.nanmean(reads_per_tn_list), ''])
                line9.append([chrom, 'Median number of reads per transposon', np.nanmedian(reads_per_tn_list), ''])
                line10.append([chrom, '25th percentile reads per transposon', np.percentile(reads_per_tn_list,25), ''])
                line11.append([chrom, '75th percentile reads per transposon', np.percentile(reads_per_tn_list,75), ''])


            elif bed_file_counter == 1:
                line1[chromosome_romannames_list.index(chrom)][-1] = len(reads_per_tn_list)
                line2[chromosome_romannames_list.index(chrom)][-1] = len(tn_insertion_position_list)/chr_length_dict.get(chrom)*100
                
                line3[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_meanfrequency
                line4[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_medianfrequency
                line5[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_25percentilefrequency
                line6[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_75percentilefrequency
                
                line7[chromosome_romannames_list.index(chrom)][-1] = max(bp_between_tn_insertions)
                
                line8[chromosome_romannames_list.index(chrom)][-1] = np.nanmean(reads_per_tn_list)
                line9[chromosome_romannames_list.index(chrom)][-1] = np.nanmedian(reads_per_tn_list)
                line10[chromosome_romannames_list.index(chrom)][-1] =  np.percentile(reads_per_tn_list,25)
                line11[chromosome_romannames_list.index(chrom)][-1] = np.percentile(reads_per_tn_list,75)
#%% DETERMINE STATISTICS FOR THE ENTIRE GENOME
#        if chromosome == None:
        bp_between_tn_insertions_genome = []
        number_tn_insertions_list = []
        reads_per_tn_genome = []
        number_tn_insertions_genome = 0
        for chrom in chrom_loop:
            number_tn_insertions_genome += len(reads_per_tn_dict.get(chrom))
            #!!!the next line includes the distance between the start of each chromosome and the first insertion and the distance between the last insertion and the end of the chromosome.
            #This might not be accurate. Please check!
            for bp_between in bp_between_tn_insertions_dict.get(chrom):
                bp_between_tn_insertions_genome.append(bp_between)
            number_tn_insertions_list.append(len(bp_between_tn_insertions_dict.get(chrom)))
            for reads_tn in reads_per_tn_dict.get(chrom):
                reads_per_tn_genome.append(reads_tn)

        if bed_file_counter == 0:
            line1.append(['Genome','Number of insertions',number_tn_insertions_genome, ''])
            line2.append(['Genome', 'Coverage percent', sum(number_tn_insertions_list)/sum(chr_length_dict.values())*100, ''])
            line3.append(['Genome', 'Mean distance between insertions', np.nanmean(bp_between_tn_insertions_genome), ''])
            line4.append(['Genome', 'Median distance between insertions', np.nanmedian(bp_between_tn_insertions_genome), ''])
            line5.append(['Genome', '25th percentile distance between insertions', np.percentile(bp_between_tn_insertions_genome,25), ''])
            line6.append(['Genome', '75th percentile distance between insertions', np.percentile(bp_between_tn_insertions_genome,75), ''])
            line7.append(['','','', ''])
            line8.append(['Genome', 'Mean number of reads per transposon', np.nanmean(reads_per_tn_genome), ''])
            line9.append(['Genome', 'Median number of reads per transposon', np.nanmedian(reads_per_tn_genome), ''])
            line10.append(['Genome', '25th percentile reads per transposon', np.percentile(reads_per_tn_genome,25), ''])
            line11.append(['Genome', '75th percentile reads per transposon', np.percentile(reads_per_tn_genome,75), ''])

        elif bed_file_counter == 1:
            line1[-1][-1] = number_tn_insertions_genome
            line2[-1][-1] = sum(number_tn_insertions_list)/sum(chr_length_dict.values())*100
            line3[-1][-1] = np.nanmean(bp_between_tn_insertions_genome)
            line4[-1][-1] = np.nanmedian(bp_between_tn_insertions_genome)
            line5[-1][-1] = np.percentile(bp_between_tn_insertions_genome,25)
            line6[-1][-1] = np.percentile(bp_between_tn_insertions_genome,75)
            line7[-1][-1] = ''
            line8[-1][-1] = np.nanmean(reads_per_tn_genome)
            line9[-1][-1] = np.nanmedian(reads_per_tn_genome)
            line10[-1][-1] = np.percentile(reads_per_tn_genome,25)
            line11[-1][-1] = np.percentile(reads_per_tn_genome,75)

        bed_file_counter += 1
#%% WRITE TO TEXT FILE
    print('Writing to text file...')

    header0 = ['chromosome','item','Dataset 1','Dataset 2']
    header  = ['          ','    ','         ','         ']

    with open(text_file,'a') as t:
        for i in range(0,len(line1)):
            table = [line1[i], line2[i], line3[i], line4[i], line5[i],
                     line6[i], line7[i], line8[i], line9[i],
                     line10[i], line11[i]]
            if i == 0:
                t.write(tabulate(table,tablefmt='github',headers=header0))
            else:
                t.write(tabulate(table,tablefmt='github',headers=header))
            t.write('\n')

#%%
if __name__ == '__main__':
    textfile_dataset_compare(bed_files=[r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed",
                                      r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_Trimmed_Aligned\Cerevisiae_WT1_Michel2017_Trimmed_Aligned.sorted.bam.bed"],
                                     text_file=r"C:\Users\gregoryvanbeek\Desktop\test_DataCompare.txt")
#    textfile_dataset_compare(bed_files=[r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed"],
#                                     text_file=r"C:\Users\gregoryvanbeek\Desktop\test_DataCompare.txt")