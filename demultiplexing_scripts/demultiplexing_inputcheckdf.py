# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 07:43:46 2020

@author: gregoryvanbeek

This is a code to test the output of demulitplexing.py
"""

#"GCCACATA"
#"GCGAGTAA"
#"GAGCTGAA"
#"GATAGACA"

#"GCTTAGGTAGTGCTGGTCGTAGTGAGCCACATAtccgtcccgcaagttaaata"
#"AGGTCAGTCACATGGTTAGGACGCAGCGAGTAAacgaaaacgaacgggataaa"
#"GCTTAGGTAGTGCTGGTCGTAGTGAGAGCTGAAtccgtcccgcaagttaaata"
#"AGGTCAGTCACATGGTTAGGACGCAGATAGACAacgaaaacgaacgggataaa"
#
#files_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
#              r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"]
#
#barcode_list = ["GCCACATA",
#                "GCGAGTAA",
#                "GAGCTGAA",
#                "GATAGACA"]

#sample_list = [1,1,2,2]

import os
import pandas as pd

def which_barcode_in_read(files_list=[], barcode_list=[], sample_list=[], saving_dataframe=False):
    '''
    Get the barcode in all reads in the files and counts how often each barcode is present in all files.
    
    files_list: list of paths to fastq files.
    barcode_list: list of all barcodes.
    sample_list: list of samples to which the barcodes belong. This must be in the same order as the barcodes.
    '''

    if len(files_list) == 0:
        print("No input files are given. Required input: files_list[file1.fq, file2.fq, ...]")
        return
    else:
        for file in files_list:
            if not os.path.isfile(file):
                print("File %s does not exist." % file)
                return

    if len(barcode_list) == 0:
        print("No barcodes are given. Required input: barcode_list[barcode1, barcode2, ...]")
        return
    if len(sample_list) == 0:
        print("No samples are given. Required input: sample_list[sample_number1, sample_number2, ...]")
        return
    elif not len(barcode_list) == len(sample_list):
        print("barcode_list and sample_list do not have the same length. Enter just as many sample numbers as there are barcode sequences.")
        return
    


    barcode_counter_dict = {}
    for barcode in barcode_list:
        barcode_counter_dict[barcode] = 0
    del (barcode)



    barcode_sample_dict = {}
    ind = 0
    for barcode in barcode_list:
        barcode_sample_dict[barcode] = sample_list[ind]
        ind += 1
    del (barcode, sample_list, ind)



    line_counter = 0
    file_counter = 1
    read_header_list = []
    read_seq_list = []
#    read_phred_list = []
    read_barcode_list = []
    read_sample_list = []
    for file in files_list:
        with open(file, 'r') as f:
            print("Analyzing %s" % os.path.basename(file))

            for line in f:
                if line.startswith('@'):
                    line_counter += 1
#                    if line_counter % 100000:
#                        print(line_counter)

                    line_header = line
                    read_header_list.append(line_header)
                    line_seq = f.readline()
                    read_seq_list.append(line_seq)
                    line_dummy = f.readline()
                    if not line_dummy[0] == '+':
                        print('WARNING: dummy line in sequence is not +.')
                        print(line_header)
#                    line_phred = f.readline()
#                    read_phred_list.append(line_phred)

                    barcode_counter = 0
                    for barcode in barcode_counter_dict:
                        barcode_counter += 1

                        if barcode_counter <= len(barcode_list) and barcode.upper() in line_seq:
                            barcode_counter_dict[barcode] += 1

                            read_barcode_list.append(barcode)
                            read_sample_list.append(barcode_sample_dict.get(barcode))
                            break
                        elif barcode_counter >= len(barcode_list):
                            read_barcode_list.append('N')
                            read_sample_list.append(0)

#                    if line_counter > 10:
#                        break
#                    else:
#                        print(line_header)
#                        print(line_seq)

        file_counter += 1

    reads_dict = {"header": read_header_list,
                  "sample": read_sample_list}

    reads_df = pd.DataFrame(reads_dict, columns = [column_name for column_name in reads_dict])
    
    reads_df.sort_values(by=['header'], ignore_index=True, inplace=True)


    pairsample_list = [0]*len(reads_df)
    unknown_reads_list = []
    for index_df, row in reads_df.iterrows():
        if row['header'].split(' ')[0] == reads_df.loc[index_df+1, 'header'].split(' ')[0]:
            if row['sample'] == reads_df.loc[index_df+1, 'sample'] and not row['sample'] == 0:
                pairsample_list[index_df] = row['sample']
                pairsample_list[index_df+1] = row['sample']
            elif not row['sample'] == reads_df.loc[index_df+1, 'sample'] or row['sample'] == 0 or reads_df.loc[index_df+1, 'sample'] == 0:
                unknown_reads_list.append(row['header'])
                unknown_reads_list.append(reads_df.loc[index_df+1, 'header'])
                reads_df.drop(reads_df.index[index_df], inplace=True)
                reads_df.drop(reads_df.index[index_df+1], inplace=True)
                del pairsample_list[index_df]
                del pairsample_list[index_df+1]

    reads_df['pairsample'] = pairsample_list

#    reads_to_sample1_list = []
#    reads_to_sample2_list = []
#    reads_sample12_list = []
#    reads_sample21_list = []
#    reads_unknown_list = []
#    for index, row in reads_df.iterrows():
#        if row['header'].split(' ')[0] == reads_df.loc[index+1, 'header'].split(' ')[0]:
#            if row['sample'] == reads_df.loc[index+1, 'sample'] and row['sample'] == 1:
#                reads_to_sample1_list.append(row['header'])
#                reads_to_sample1_list.append(reads_df.loc[index+1, 'header'])
#            elif row['sample'] == reads_df.loc[index+1, 'sample'] and row['sample'] == 2:
#                reads_to_sample2_list.append(row['header'])
#                reads_to_sample2_list.append(reads_df.loc[index+1, 'header'])
#            elif not row['sample'] == reads_df.loc[index+1, 'sample'] and row['sample'] == 1 and reads_df.loc[index+1, 'header'] == 2:
#                reads_sample12_list.append(row['header'])
#                reads_sample12_list.append(reads_df.loc[index+1, 'header'])
#            elif not row['sample'] == reads_df.loc[index+1, 'sample'] and row['sample'] == 2 and reads_df.loc[index+1, 'header'] == 1:
#                reads_sample21_list.append(row['header'])
#                reads_sample21_list.append(reads_df.loc[index+1, 'header'])
#            elif row['sample'] == 0 or reads_df.loc[index+1, 'sample'] == 0:
#                reads_unknown_list.append(row['header'])
#                reads_unknown_list.append(reads_df.loc[index+1, 'header'])


#    pair_sample = [0]*len(reads_df)
#    for index, row in reads_df.iterrows():
##        print('Current index: %i with sample number %i' % (index, row['sample']))
#        header_split = row['header'].split(' ') #split header to find pair number (first number after space in header)
#        if header_split[1].startswith('1'): #only analyze pair number 1. next line searches for same header for pair number 2.
#            readpair_sample = reads_df.loc[reads_df['header'] == header_split[0] + ' ' + '2' + header_split[1][1:], 'sample'] #search for the same header, but with pair number 2
#            if len(readpair_sample) > 0: #check whether pair is found
##                print('sample read pair = %i ' % readpair_sample.iloc[0])
#                if row['sample'] == readpair_sample.iloc[0]: #check whether sample numbers for pair number 1 and pair number 2 are the same.
##                    print('samples are the same.')
#                    pair_sample[index] = row['sample'] #put sample number for pair number 1
#                    pair_sample[readpair_sample.index.values[0]] = readpair_sample.iloc[0] #put sample number for pair number 2
#                else: #for cases where sample numbers are not the same.
##                    print('samples are different.')
#                    if row['sample'] < readpair_sample.iloc[0]: #if pair number 1 belongs to sample 1 and pair 2 to sample 2
#                        pair_sample[index] = -1
#                        pair_sample[readpair_sample.index.values[0]] = -1
#                    elif row['sample'] > readpair_sample.iloc[0]: #if pair number 1 belongs to sample 2 and pair 2 to sample 1
#                        pair_sample[index] = -2
#                        pair_sample[readpair_sample.index.values[0]] = -2
#
#
#    reads_df['pair_belong_to_sample'] = pair_sample


    if saving_dataframe == True:
        reads_df.to_csv(os.path.join(os.path.dirname(files_list[0]), 'reads_dataframe.csv'), index=False)
    
    return(reads_df)

#    unmatched_reads = abs(sum([hits for key, hits in barcode_counter_dict.items()]) - line_counter)
#    barcode_counter_dict['unmatched_reads'] = unmatched_reads
#
#    print('')
#    print("Total number of reads found: %i" % line_counter)
#    print(barcode_counter_dict)
#
#    print('First reads:')
#    i = 0
#    for k,v in read_barcode_dict.items():
#        if i <= 5:
#            print(k)
#            print(v)
#            print('')
#        else:
#            break
#        i += 1




#barcode_sample_dict = {"GCCACATA": 1,
#                       "GCGAGTAA": 1,
#                       "GAGCTGAA": 2,
#                       "GATAGACA": 2}
#
#for read, barcode in read_barcode_dict.items():
#    if read.split(' ')[1][0] == '1':
#        barcode1 = barcode
#        barcode2 = read_barcode_dict.get(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
#        if not barcode_sample_dict.get(barcode1) == barcode_sample_dict.get(barcode2) and not barcode2 == None:
#            print(read)
#            print(read.split(' ')[0] + ' 2' + read.split(' ')[1][1:])
#            print(barcode1)
#            print(barcode2)
#            break
#            pass
        


if __name__ == '__main__':
    reads_df = which_barcode_in_read(files_list = [r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
                                        r"C:\Users\gregoryvanbeek\Documents\Data_Sets\20201104_Enzo\raw_data_20201104\D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"],
                          barcode_list = ["GCCACATA",
                                          "GCGAGTAA",
                                          "GAGCTGAA",
                                          "GATAGACA"],
                          sample_list = [1,1,2,2],
                          saving_dataframe=True)
