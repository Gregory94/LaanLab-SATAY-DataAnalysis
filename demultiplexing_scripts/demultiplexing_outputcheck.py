# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 13:46:31 2020

@author: gregoryvanbeek
"""


import os, sys



def check_fastq_pairedinterleaved(inputfile='', barcodes=[]):

    if not os.path.isfile(inputfile):
        print('WARNING: input file not found.')
        sys.exit()
    if barcodes == []:
        print('WARNING: barcodes list is empty')
        sys.exit()



    if not barcodes == []:
        print("Assessing barcodes: ", barcodes)



    total_line_counter = 0
    header_counter = 0
    with open(inputfile,'r') as f:
        for line in f:
            total_line_counter += 1
            if line.startswith('@'):
                header_counter += 1

    print("Number of lines = %i" % total_line_counter)
    print("Number of reads = %i" % header_counter)
    lines_divisible_by_four = (total_line_counter / 4).is_integer()
    if lines_divisible_by_four == False:
        print("ERROR: Number of lines in fastq file is not a multiple of 4. Format of fastq file may be corrupted.")
        sys.exit()
    
    lines_divisible_by_headercounter = total_line_counter / header_counter
    if not int(lines_divisible_by_headercounter) == 4:
        print("ERROR: Number of header in fastq file does not match the number of lines in the file. Format of fastq file may be corrupted.")
        sys.exit()




    pair_counter = 0
    with open(inputfile,'r') as f:
        for line in f:
            if line.startswith('@'):
                read_dict = {}
                read_dict['header'] = line.strip('\n')


                #Check if two subsequent headers have the same header id
                if pair_counter == 0:
                    header1 = line.split(' ')[0]
                    pair_counter += 1
                elif pair_counter == 1:
                    header2 = line.split(' ')[0]
                    if not header1 == header2:
                        print('Header %s is not the same as the previous header %s' % (header2, header1))
                    pair_counter = 0

                read_dict['sequence'] = next(f).strip('\n')
                read_dict['dummy'] = next(f).strip('\n')
                read_dict['phred'] = next(f).strip('\n')


                #Check if there are other characters than A,C,T,G or N in the sequence line.
                valid_char = "ATCGN"
                if not all(c in valid_char for c in read_dict.get('sequence')):
                    print("WARNING: Sequence contains invalid characters at header %s " % read_dict.get('header'))
#                elif 'N' in read_dict.get('sequence'):
#                    print("WARNING: Sequence contains N at header %s " % read_dict.get('header'))


                #Check if sequence and phred (quality score) have the same number of characters
                if not len(read_dict.get('sequence')) == len(read_dict.get('phred')):
                    print("WARNING: sequence and phred lines are not the same length at header %s " % read_dict.get('header'))


                #Check if barcodes are present in sequence
                if not barcodes == []:
                    for barcode in barcodes:
                        if barcode.upper() in read_dict.get('sequence'):
                            barcode_in_line = True
                            break
                    if not barcode_in_line == True:
                        print('No barcode found in read %s' % header1)
                    else:
                        barcode_in_line = False




if __name__ == '__main__':
    check_fastq_pairedinterleaved(inputfile=r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/wt1_enzo_dataset_demultiplexed_interleaved_sample1/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample1interleavedsorted_pairs.fq",
                              barcodes=["GCCACATA", "GCGAGTAA"])
#    check_fastq_pairedinterleaved(inputfile=r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/wt1_enzo_dataset_demultiplexed_interleaved/wt1_enzo_dataset_demultiplexed_interleaved_sample2/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_pairs.fq",
#                              barcodes=["GAGCTGAA", "GATAGACA"])
#    check_fastq_pairedinterleaved(inputfile=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\wt1_dataset_enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_sample2interleavedsorted_pairs_400lines.fq",
#                                  barcodes=["GAGCTGAA", "GATAGACA"])
