# -*- coding: utf-8 -*-
"""
"""

#%%
def chromosome_name_bedfile(bed_file=None):
    '''This function returns some properties of the chromosomes in a bed file.
    Input can be either of two options:
        The full path to a bed file after which the program opens the bed file, or
        a list of the lines in the bed file. The latter requires to read the bed file before calling this function and input all lines in the bed file as a list. The function than does not open the bed file again.
    Returns three dictionaries (in this order):
        The first indicates the names of the chromosomes as used in the bed file (keys are roman numerals 1 to 16 and the values are the names used in the bed file).
        The second is the start line in the bed file of each chromosome (keys are the roman numerals of the chromosome names and the values are the start lines in the bed file of the chromosome).
        The third is the end line in the bed file of each chromosome (keys are the roman numerals of the chromosome names and the values are the start lines in the bed file of the chromosome)
    '''
    
    if bed_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        bed_file = os.path.join(file_dirname,'..','..','satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed1.bam.bed')



    
    if type(bed_file) is str:
        with open(bed_file) as f:
            lines = f.readlines()
    elif type(bed_file) is list:
        lines = bed_file




    num_arabic = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    num_roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
    roman_to_arabic_dict = {}
    index_counter = 0
    for rom in num_roman:
        roman_to_arabic_dict[rom] = num_arabic[index_counter]
        index_counter += 1

    chromosome_romannames_list = []
    for roman in roman_to_arabic_dict:
        chromosome_romannames_list.append(roman)




    chrom_names_dict = {}
    chrom_start_line_dict = {}
    chrom_end_line_dict = {}

    chrom_name_in_bed = ''
    chr_counter = 0
    line_counter = 0
    stop_loop = False
    while stop_loop is False:
        line = lines[line_counter]
        chrom_name_current = line.split(' ')[0].replace('chr','')
        if not chrom_name_current.startswith('track') and not chrom_name_current.startswith('M'): #SKIP HEADER AND MITOCHRONDRIAL CHROMOSOMES
            if chrom_name_current != chrom_name_in_bed:
                chrom_names_dict[chromosome_romannames_list[chr_counter]] = chrom_name_current
                chrom_name_in_bed = chrom_name_current
#                print('Chromosome ',chromosome_romannames_list[chr_counter], 'is ',chrom_name_current)
                
                chrom_start_line_dict[chromosome_romannames_list[chr_counter]] = line_counter #GET START INDEX IN THE BED FILE OF THE CURENT CHROMOSOME
                if chr_counter != 0:
                    chrom_end_line_dict[chromosome_romannames_list[chr_counter-1]] = line_counter-1 #GET THE END INDEX IN THE BED OF THE PREVIOUS CHROMOSOME (SKIP FOR THE FIRST CHROMOSOME)

                chr_counter += 1

        elif chrom_name_current.startswith('M'):
            chrom_end_line_dict[chromosome_romannames_list[-1]] = line_counter-1 #GET THE END INDEX IN THE BED FILE FOR THE FINAL CHROMOSOME
            stop_loop = True
                
        line_counter += 1
        
    return(chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict)

#%%
def chromosome_name_wigfile(wig_file=None):
    if wig_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        wig_file = os.path.join(file_dirname,'..','..','satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed1.bam.wig')






    if type(wig_file) is str:
        with open(wig_file) as f:
            lines = f.readlines()
    elif type(wig_file) is list:
        lines = wig_file




    num_arabic = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    num_roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
    roman_to_arabic_dict = {}
    index_counter = 0
    for rom in num_roman:
        roman_to_arabic_dict[rom] = num_arabic[index_counter]
        index_counter += 1

    chromosome_romannames_list = []
    for roman in roman_to_arabic_dict:
        chromosome_romannames_list.append(roman)





    chrom_names_dict = {}
    chrom_start_line_dict = {}
    chrom_end_line_dict = {}

    chr_counter = 0
    line_counter = 0
    for line in lines:
        line.strip('\n')
        chrom_line = 'variableStep'
        line_split = line.split(' ')
        if line_split[0] == chrom_line:
            chromosome_name = line_split[1].replace('chrom=chr','').strip('\n')
            chrom_names_dict[chromosome_romannames_list[chr_counter]] = chromosome_name
            print('Chromosome ',chromosome_romannames_list[chr_counter], 'is ',chromosome_name)
            
            chrom_start_line_dict[chromosome_romannames_list[chr_counter]] = line_counter+2 #GET START INDEX IN THE BED FILE OF THE CURENT CHROMOSOME
            if chr_counter != 0:
                chrom_end_line_dict[chromosome_romannames_list[chr_counter-1]] = line_counter #GET THE END INDEX IN THE BED OF THE PREVIOUS CHROMOSOME (SKIP FOR THE FIRST CHROMOSOME)

            chr_counter += 1
        line_counter += 1

    chrom_end_line_dict[chromosome_romannames_list[chr_counter-1]] = len(lines)


    return(chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict)

#%%
if __name__ == '__main__':
    chromosome_name_wigfile()