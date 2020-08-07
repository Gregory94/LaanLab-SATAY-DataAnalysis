# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 13:41:44 2020

@author: gregoryvanbeek



This python code creates a .txt file that includes all essential genes that are annotated as essential genes in Saccharomyces Cerevisiae.
It gets this information from multiple .txt files and combines all information in one file.

It checks whether a gene is already present in the file (either with the same name or an alias).
It stores the gene names using the oln naming convention.
"""



import os, sys

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from gene_names import gene_aliases


#%%
def create_essentialgenes_list(inputfiles_list = None):
    '''
    This function requires as input a list of paths to files containing essential genes.
    Multiple files can be present in this list.
    The input files have to have the following layout:
        - Three header lines (can be empty or containing any text)
        - Each new lines should contain one gene name in either oln or designation naming convention.
    
    This function is dependable on the following custom made modules:
        - gene_names.py (requires the file Yeast_Protein_Names.txt)
    
    The output will be a text file containing all uniquely found genes in all input files given.
    The file will be stored at the same location of the first file of the input list with the name 'Cerevisiae_AllEssentialGenes_List.txt'.
    In this file each line contains one gene and it has a single header line containing all the filenames that were used to create this file.
    
    '''
    
    if inputfiles_list == None:
        ValueError('Input list containing one or more paths is missing.')
    else:
        files = inputfiles_list


    path = os.path.dirname(files[0])
    filename_list = []
    for file in files:
        filename_list.append(os.path.basename(file))

    del (inputfiles_list, file)
    
    
#%%
    all_genes_list = []
    for file in files:#ASSUMES THREE HEADER LINES
        filename = os.path.basename(file)
        with open(file) as f:
            lines = f.readlines()
            print('Number of genes found in %s: %i' % (filename, (len(lines)-3)))

            for line in lines[3:]:
                all_genes_list.append(line.strip('\n'))


    del (file, f, lines, line)


#%%
    gene_aliases_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Yeast_Protein_Names.txt")[0]


#%%
    all_genes_oln_list = []
    for gene in all_genes_list:
        if gene in gene_aliases_dict:
            all_genes_oln_list.append(gene)
        else:
            for key, val in gene_aliases_dict.items():
                if gene in val:
                    all_genes_oln_list.append(key)
                    break


    del (gene, all_genes_list, key, val, gene_aliases_dict)


#%%
    unique_genes_list = list(set(all_genes_oln_list))
    unique_genes_list.sort()

    print('Number of unique essential genes found : %i' % len(unique_genes_list))


    del (all_genes_oln_list)


#%%
    save_filename = r'Cerevisiae_AllEssentialGenes_List.txt'
    save_file = os.path.join(path, save_filename)

    print('Creating text file with all unique genes at %s' % save_file)


    with open(save_file, 'w') as f:
        f.write('All essential genes found in lists:' + str(filename_list) + '\n')
        for gene in unique_genes_list:
            f.write(gene + '\n')

    del (gene)

#%%
if __name__ == '__main__':
    create_essentialgenes_list([r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_EssentialGenes_List_1.txt",
                 r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_EssentialGenes_List_2.txt"])