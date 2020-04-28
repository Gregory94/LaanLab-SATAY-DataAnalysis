# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 09:25:54 2020

@author: gregoryvanbeek
"""

import sys
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
sys.path.insert(1,r'C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\python_modules')
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic

#%%
def chromosome_insertion_periodicity(chromosome=None,bed_file=None):
    '''
    '''

#%% USED FILES
    gff_file = r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3"

#%% GET CHROMOSOME START AND END POSTIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)

#%% GET ROMAN ARABIC NUMERALS
    roman_to_arabic_dict = chromosomename_roman_to_arabic()[1]
    chromosome_romannames_list = []
    for roman in roman_to_arabic_dict:
        chromosome_romannames_list.append(roman)

#%% OPEN BED FILE
    with open(bed_file) as f:
        lines = f.readlines()

#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
    chrom_names_dict = {}
    chrom_start_index_dict = {}
    chrom_end_index_dict = {}
    chrom_name = ''
    chr_counter = 0
    line_counter = 0
    stop_loop = False
    while stop_loop is False:
        line = lines[line_counter]
        chrom_name_current = line.split(' ')[0].replace('chr','')
        if not chrom_name_current.startswith('track'): #SKIP HEADER
            if not chrom_name_current.startswith('M'): #SKIP MITOCHRONDRIAL CHROMOSOMES
                if chrom_name_current != chrom_name:
                    chrom_names_dict[chromosome_romannames_list[chr_counter]] = chrom_name_current
                    chrom_name = chrom_name_current
#                    print('Chromosome ',chromosome_romannames_list[chr_counter], 'is ',chrom_name_current)
                    
                    chrom_start_index_dict[chromosome_romannames_list[chr_counter]] = line_counter #GET START INDEX IN THE BED FILE OF THE CURENT CHROMOSOME
                    if chr_counter != 0:
                        chrom_end_index_dict[chromosome_romannames_list[chr_counter-1]] = line_counter-1 #GET THE END INDEX IN THE BED OF THE PREVIOUS CHROMOSOME (SKIP FOR THE FIRST CHROMOSOME)

                    chr_counter += 1

            elif chrom_name_current.startswith('M'):
                chrom_end_index_dict[chromosome_romannames_list[-1]] = line_counter-1 #GET THE END INDEX IN THE BED FILE FOR THE FINAL CHROMOSOME
                stop_loop = True
                
        line_counter += 1

#%% DETERMINE STATISTICS FOR INDIVIDUAL CHROMOSOMES AND PUT THE RESULTS IN A DICT
    if chromosome != None:
        chrom_loop = {}
        chrom_loop[chromosome] = chrom_names_dict.get(chromosome)
    else:
        chrom_loop = chrom_names_dict

    bp_between_tn_insertions_dict = {}
    for chrom in chrom_loop:
        tn_insertion_position_list = []
        for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
            line = line.strip('\n').split()
            tn_insertion_position_list.append(int(line[1]))
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(tn_insertion_position_list[:-1], tn_insertion_position_list[1:])]
        bp_between_tn_insertions.insert(0,tn_insertion_position_list[0]) #ADD START OF GENE (bp=0)
        bp_between_tn_insertions.append(chr_length_dict.get(chrom) - tn_insertion_position_list[-1]) #ADD END OF GENE (bp=INDEX LAST TN - GENE LENGTH)
        bp_between_tn_insertions_dict[chrom] = bp_between_tn_insertions

        print('')
        print('For chromosome ',chrom,':')
        print('Coverage is %.2f percent' % (len(tn_insertion_position_list)/chr_length_dict.get(chrom)*100))
        tn_insertion_meanfrequency = np.nanmean(bp_between_tn_insertions)
        tn_insertion_medianfrequency = np.nanmedian(bp_between_tn_insertions)
        print('Mean transposon insertion periodicity is once every %.2f bp' % tn_insertion_meanfrequency)
        print('Median transposon insertion periodicity is once every %.2f bp' % tn_insertion_medianfrequency)
        print('')

#%% APPLY AUTOCORRELATION FOR CHECKING THE PERIODICITY

    #MAKE LIST OF ALL INSERTION LOCATIONS AND HOW MANY INSERTIONS EACH LOCATION HAS
#    for chrom in chrom_loop:
#        number_of_insertions = [0]*chr_length_dict.get(chrom)
#        for ins in tn_insertion_position_list:
#            number_of_insertions[ins] += 1
#    
#    norm = number_of_insertions - np.mean(number_of_insertions)
#    n = norm.size
#    corr = np.correlate(norm, norm, mode='same')
#    autocorr = corr[n//2 + 1:] / (np.var(number_of_insertions) * np.arange(n-1, n//2, -1))
#    lag = np.abs(autocorr).argmax() + 1
#    print(lag)
#    r = autocorr[lag-1]
#    print(r)
    
#%% DETERMINE STATISTICS FOR THE ENTIRE GENOME
    if chromosome == None:
        bp_between_tn_insertions_genome = []
        number_tn_insertions_list = []
        for chrom in chrom_loop:
            #the next line includes the distance between the start of each chromosome and the first insertion and the distance between the last insertion and the end of the chromosome.
            #This might not be accurate. Please check!
            for bp_between in bp_between_tn_insertions_dict.get(chrom):
                bp_between_tn_insertions_genome.append(bp_between)
            number_tn_insertions_list.append(len(bp_between_tn_insertions_dict.get(chrom)))
#            number_tn_insertions_list.append(sum(x > 0 for x in alltransposoncounts_dict.get(chrom)))
            
        print('')
        print('For the entire genome:')
        print('Coverage is %.2f percent' % (sum(number_tn_insertions_list)/sum(chr_length_dict.values())*100))
        print('Mean transposon insertion periodicity for the entire genome is %.2f' % np.nanmean(bp_between_tn_insertions_genome))
        print('Median transposon insertion periodicity for the entire genome is %.2f' % np.nanmedian(bp_between_tn_insertions_genome))

#%% DETERMINE THE DISTRIBUTION OF THE NUMBER OF BP BETWEEN SUBSEQUENT TRANSPOSON INSERTIONS

    if chromosome != None:
        df = pd.DataFrame(data=bp_between_tn_insertions_dict)
        df.columns = [chromosome]
        df_melt = df.melt(var_name='chromosomes',value_name='bp between insertions')
    elif chromosome == None:
        df = pd.DataFrame(data=bp_between_tn_insertions_dict.get('I'))
        for chrom in chrom_loop:
            if chrom != 'I':
                df_temp = pd.DataFrame(data=bp_between_tn_insertions_dict.get(chrom))
                df = pd.concat([df,df_temp], axis=1)
        df.columns = chromosome_romannames_list
        df_melt = df.melt(var_name='chromosomes',value_name='bp between insertions')
    
    sb.violinplot(x='chromosomes',y='bp between insertions',data=df_melt,inner='box',gridsize=1000, cut=0)
#    v.set_yscale('log')

#%%
    return(tn_insertion_meanfrequency,bp_between_tn_insertions_dict)

#%%
if __name__ == '__main__':
    chromosome_insertion_periodicity(chromosome='I',bed_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed")