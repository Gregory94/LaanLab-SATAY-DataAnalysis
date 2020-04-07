# -*- coding: utf-8 -*-
"""
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(1,r'python_modules')
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic


def transposon_profile(chrom='I',bar_width=1000,binsize=1000,wig_file = None):
    '''This function creates a bar plot along a specified chromosome for the number of transposons.
    The height of each bar represents the number of transposons at the position indicated on the x-axis.
    The input is as follows: which chromosome (indicated by roman numeral), bar_width, binsize, wig_file.
    The bar_width determines how many basepairs are put in one bin. little basepairs per bin may be slow. Too many basepairs in one bin and possible low tranposon areas might be obscured.
    The binsize  is only for visualization, choose a value such that the plot has a high constrast in transposon rich areas and transposon poor areas.
    The wig_file is one of the files created by the Matlab code from the kornmann-lab.
    '''
    #GET CHROMOSOME LENGTHS AND POSITIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3")
    
    
    #CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
    arabic_to_roman_dict, roman_to_arabic_dict = chromosomename_roman_to_arabic()    
    chromosomenames_list = []
    for roman in roman_to_arabic_dict:
        chromosomenames_list.append(roman)
    
    
#    for chrom in chromosomenames_list:
    chrom_index = chromosomenames_list.index(chrom)
    
#    allinsertionsites_list = list(range(0,chr_length_dict.get(chrom))) #CREATE LIST OF ALL POSIBLE INSERTION SITES IN THE CURRENT CHROMOSOME
    alltransposoncounts_list = np.zeros(chr_length_dict.get(chrom)) #FOR EACH INSERTION SITE LIST THE NUMBER OF TRANSPOSON INSERTION. BY DEFAULT THIS 0 AND IS LATER UPDATED IF AN INSERTION SITE IS PRESENT IN THE WIG FILE
    print('Chromosome length: ',chr_length_dict.get(chrom))
    
    with open(wig_file) as f:
        lines = f.readlines()
    
    #GET ALL LINES WITH THE TRANSPOSONS FOR THE CURRENT CHROMOSOME
    line_counter = 0
    for line in lines:
        line = line.strip('\n')
        if line.endswith('chrom=chr'+chromosomenames_list[chrom_index]):
            wigfile_start_index = line_counter + 1
        elif chromosomenames_list[chrom_index] == chromosomenames_list[-1]: #CHECK IF THE LAST CHROMOSOME IS REACHED, SINCE THEN THE NEXT CHROMOSOME DOES NOT NEED TO BE SEARCHED AS THIS WON'T EXISTS
            wigfile_end_index = len(lines)-1 #GET INDEX LAST ELEMENT
        elif line.endswith('chrom=chr'+chromosomenames_list[chrom_index+1]):
            wigfile_end_index = line_counter
        line_counter += 1
    


    #GET ALL TRANSPOSON COUNTS FOR THE CURRENT CHROMOSOME
    for line in lines[wigfile_start_index:wigfile_end_index]:
        line = line.strip(' \n').split()
        alltransposoncounts_list[int(line[0])] = int(line[1])
    
    
    
    
    #THE LIST WITH ALL THE TRANPOSONS FOR THE CURRENT CHROMOSOME IS TYPICALLY REALLY LARGE.
    #TO COMPRESS THIS LIST, THE BASEPAIR POSITIONS ARE GROUPED IN GROUPS WITH SIZE DEFINED BY 'BAR_WIDTH'
    #IN EACH GROUP THE NUMBER OF TRANSPOSONS ARE SUMMED UP.
    #THIS IS DONE TO SPEED UP THE SCRIPT AS PLOTTING ALL VALUES IS SLOW
#    bar_width = 1000
    alltransposoncounts_binnedlist = []
    val_counter = 0
    sum_values = 0
    if bar_width == 1:
        alltransposoncounts_binnedlist = alltransposoncounts_list
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/float(bar_width)))
    else:
        for n in range(len(alltransposoncounts_list)):
            if val_counter % bar_width != 0:
                sum_values += alltransposoncounts_list[n]
            elif val_counter % bar_width == 0:
                alltransposoncounts_binnedlist.append(sum_values)
                sum_values = 0
            val_counter += 1
            
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/bar_width)+1)
    
    
    print('Plotting chromosome ', chrom, '...')
    fig,ax = plt.subplots()
    ax.bar(allinsertionsites_list,alltransposoncounts_binnedlist,width=binsize)
    ax.set_yscale('log')
    plt.show()
    

if __name__ == '__main__':
    transposon_profile('XVI',1000,1000,r"E-MTAB-4885.WT1.bam.wig")