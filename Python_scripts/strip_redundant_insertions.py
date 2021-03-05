# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:39:53 2021

@author: gregoryvanbeek
"""

import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position
from chromosome_names_in_files import chromosome_name_bedfile




chr_length_dict = chromosome_position()[0]

filepath = r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dnrp1.bed"

filepath_splitext = os.path.splitext(filepath)
exten = filepath_splitext[1]




if exten == ".bed":
    print("Bed file loaded %s" % filepath)
    
    chrom_start_line_dict, chrom_end_line_dict = chromosome_name_bedfile(filepath)[1:3]
    
    with open(filepath, 'r') as f:
        lines = f.readlines()


    with open(filepath_splitext[0]+"_clean.bed", "w") as w:
        w.write(lines[0])
        for chrom in chr_length_dict:
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
else:
    print("Extension not recognized")