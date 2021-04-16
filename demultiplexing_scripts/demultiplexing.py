# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a python script for demultiplexing reads from transposon sequencing data.
"""

import os


#DEFINE INPUT FILES (UNZIPPED FASTQ FORMAT)
inputfile_list = [r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1.fq",
                  r"/home/laanlab/Documents/satay/datasets/wt1_enzo_dataset/D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_2.fq"]


#DEFINE BARCODES. KEEP THE SAME NAMING FORMAT FOR THE SAMPLES (I.E. sample#_...)
barcode_dict = {"GCCACATA": 'sample1_forward',
                "GCGAGTAA": 'sample1_reverse',
                "GAGCTGAA": 'sample2_forward',
                "GATAGACA": 'sample2_reverse'}



#DEFINE OUTPUT FILESNAMES
print("Creating files for storing output data...")
print('')
outputfile_dict = {}
for barcode, samplename in barcode_dict.items():
    outputfile_dict[samplename] = os.path.splitext(inputfile_list[0])[0] + '_' + samplename + ".fq"
outputfile_dict['unmatched'] = os.path.splitext(inputfile_list[0])[0] + "_unmatched.fq"
print(outputfile_dict)
print('')



#CREATE OUTPUT FILES
for samplename, outputfile in outputfile_dict.items():
    if not os.path.exists(outputfile):
        open(outputfile,'x')



#OPEN INPUT FILES
f1 = open(inputfile_list[0],'r')
f2 = open(inputfile_list[1],'r')


#CHECK ALL READS IN INPUT FILES FOR BARCODES AND STORE THOSE IN THE DEDICATED OUTPUT FILES.
#temp=0
sample_counter = 0
unmatchedsamples_counter = 0
nobarcode_counter = 0
for line1 in f1:
#    if temp < 10:
    line2 = f2.readline()
    if line1.startswith('@') and line2.startswith('@') and line1.split(' ')[0] == line2.split(' ')[0]:
        line1_header = line1
        line2_header = line2
        line1_seqnc = f1.readline()
        line2_seqnc = f2.readline()
        line1_dummy = f1.readline()
        line2_dummy = f2.readline()
        line1_phred = f1.readline()
        line2_phred = f2.readline()
#            print(line1_header)
#            print(line1_seqnc)
#            print(line1_dummy)
#            print(line1_phred)
#            print(line2_header)
#            print(line2_seqnc)
#            print(line2_dummy)
#            print(line2_phred)


        bc1 = [barcode for barcode in barcode_dict if barcode in line1_seqnc]
#            if not bc1 == []:
#                print(bc1[0], barcode_dict.get(bc1[0]))
#            else:
#                print(bc1)
        bc2 = [barcode for barcode in barcode_dict if barcode in line2_seqnc]
#            if not bc2 == []:
#                print(bc2[0], barcode_dict.get(bc2[0]))
#            else:
#                print(bc2)

        
        if not bc1 == [] and not bc2 == []:
            if barcode_dict.get(bc1[0]).split('_')[0] == barcode_dict.get(bc2[0]).split('_')[0]:
                with open(outputfile_dict.get(barcode_dict.get(bc1[0])),'a') as f_out:
                    f_out.write(line1_header + line1_seqnc + line1_dummy + line1_phred)
                with open(outputfile_dict.get(barcode_dict.get(bc2[0])),'a') as f_out:
                    f_out.write(line2_header + line2_seqnc + line2_dummy + line2_phred)
                sample_counter += 1

            else: #READ PAIR DO NOT BELONG TO SAME SAMPLE
                with open(outputfile_dict.get("unmatched"),'a') as f_out:
                    f_out.write(line1_header + line1_seqnc + line1_dummy + line1_phred)
                    f_out.write(line2_header + line2_seqnc + line2_dummy + line2_phred)
                unmatchedsamples_counter += 1
        else: #NO BARCODE FOUND
            with open(outputfile_dict.get("unmatched"),'a') as f_out:
                f_out.write(line1_header + line1_seqnc + line1_dummy + line1_phred)
                f_out.write(line2_header + line2_seqnc + line2_dummy + line2_phred)
            nobarcode_counter += 1

#            temp+=1
#    elif temp >= 10:
#        break

f1.close()
f2.close()

print("Number of demulitplexed reads: %i" % sample_counter)
print("Number of reads where the forward and reverse barcode did not match the same sample: %i" % unmatchedsamples_counter)
print("Number of reads in which no barcode was found: %i" % nobarcode_counter)


#bash command for removing generated files:  rm D18524C717111_BDDP200001534-1A_HJVN5DSXY_L1_1_*
