# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:39:53 2021

@author: gregoryvanbeek

This script if for removing transposon insertions in .bed and .wig files that were mapped outside the chromosomes.
"""

import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position
from chromosome_names_in_files import chromosome_name_bedfile, chromosome_name_wigfile




def strip_redundant_ins(filepath=None):
    '''
    This code reads a .bed or .wig file and removed any insertions that were mapped outside a chromosome.
    For this, it creates a new file with the same name as the inputfile with the extension _clean.bed or _clean.wig, respectively.
    This is saved at the same location as the input file.
    In this _clean file the redundant insertions that were mapped outside the chromosome are removed.
    The lengths of the chromosomes are determined the python function 'chromosome_position' which is part of the python module 'chromosome_and_gene_positions.py'.
    This module gets the lengths of the chromosomes from a .gff file downloaded from SGD (https://www.yeastgenome.org/).
    '''

    if filepath == None:
        sys.exit()
    else:
        assert os.path.isfile(filepath), 'File not found: %s' % filepath

    chr_length_dict = chromosome_position()[0]

    filepath_splitext = os.path.splitext(filepath)
    exten = filepath_splitext[1]



    num_roman = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']



    if exten == ".bed":
        print("Bed file loaded %s" % filepath)

        chrom_start_line_dict, chrom_end_line_dict = chromosome_name_bedfile(filepath)[1:3]

        with open(filepath, "r") as f:
            lines = f.readlines()


        with open(filepath_splitext[0]+"_clean.bed", "w") as w:
            w.write(lines[0])
            for chrom in num_roman:
                print("evaluating chromosome %s" % chrom)

                for line in lines[chrom_start_line_dict.get(chrom): chrom_end_line_dict.get(chrom)+1]:
                    line_list = " ".join(line.strip("\n").split()).split(" ")
                    if int(line_list[1]) > chr_length_dict.get(chrom) or int(line_list[1]) < 0:
                        print("Line removed: %s" % line)
                    else:
                        w.write(line)

            for line in lines[chrom_end_line_dict.get("XVI")+1:]:
                w.write(line)





    elif exten == ".wig":
        print("Wig file loaded %s" % filepath)

        chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(filepath)[1:3]

        with open(filepath, 'r') as f:
            lines = f.readlines()

        with open(filepath_splitext[0]+"_clean.wig", "w") as w:
            w.write(lines[0])

            for chrom in num_roman:
                print("evaluating chromosome %s" % chrom)

                w.write(lines[chrom_start_line_dict.get(chrom) - 1])
                for line in lines[chrom_start_line_dict.get(chrom): chrom_end_line_dict.get(chrom)]: #no '+1' in for loop, this is only for bed file
                    line_list = " ".join(line.strip("\n").split()).split(" ")
                    if int(line_list[0]) > chr_length_dict.get(chrom) or int(line_list[0]) < 0:
                        print("Line removed: %s" % line)
                    else:
                        w.write(line)

            for line in lines[chrom_end_line_dict.get("XVI"):]:
                w.write(line)



    else:
        print("Extension not recognized")



#%%
if __name__ == '__main__':
    strip_redundant_ins(filepath = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dnrp1.wig")
    # strip_redundant_ins(filepath = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dnrp1.bed")

