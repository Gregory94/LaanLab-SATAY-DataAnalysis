# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:17:11 2021

@author: gregoryvanbeek
"""

import os, sys
import pandas as pd
import re

def dataframe_from_pergenefile(datafile, verbose=True):

    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
    from essential_genes_names import list_known_essentials #import essential_genes_names from python modules directory


# read file
    assert os.path.isfile(datafile), 'File not found at: %s' % datafile

    with open(datafile) as f:
        lines = f.readlines()[1:] #skip header

    genenames_list = [None]*len(lines)
    tnpergene_list = [None]*len(lines)
    readpergene_list = [None]*len(lines) 

    line_counter = 0
    for line in lines:
        l = re.split(' |\t', line.strip('\n'))

        genenames_list[line_counter] = l[0]
        tnpergene_list[line_counter] = int(l[1])
        readpergene_list[line_counter] = int(l[2])

        line_counter += 1

    del (line, l, line_counter, datafile)


# determine number of reads per insertion per gene
    readperinspergene_list = [None]*len(lines)
    for i in range(len(tnpergene_list)):
        if not tnpergene_list[i] == 0:
            readperinspergene_list[i] = readpergene_list[i] / tnpergene_list[i]
        else:
            readperinspergene_list[i] = 0

    del (i)


# determine essential genes
    known_essential_gene_list = list_known_essentials(input_files=[os.path.join(file_dirname,'..','Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                                                                   os.path.join(file_dirname,'..','Data_Files','Cerevisiae_EssentialGenes_List_2.txt')],
                                                      verbose=verbose)

    geneessentiality_list = [None]*len(lines)
    for i in range(len(genenames_list)):
        if genenames_list[i] in known_essential_gene_list:
            geneessentiality_list[i] = True
        else:
            geneessentiality_list[i] = False

    del (lines, file_dirname, known_essential_gene_list, i)

# create dataframe
    read_gene_dict = {"gene_names": genenames_list,
                      "gene_essentiality": geneessentiality_list,
                      "tn_per_gene": tnpergene_list,
                      "read_per_gene": readpergene_list,
                      "Nreadsperinsrt": readperinspergene_list}

    read_gene_df = pd.DataFrame(read_gene_dict, columns = [column_name for column_name in read_gene_dict])

    del (read_gene_dict, genenames_list, geneessentiality_list, tnpergene_list, readpergene_list, readperinspergene_list)
    
    
    return(read_gene_df)