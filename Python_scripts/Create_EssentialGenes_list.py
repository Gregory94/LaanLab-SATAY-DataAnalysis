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


saving = True


#%%
path = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files"

filename_list = ['Cervisiae_EssentialGenes_List_1.txt', 'Cervisiae_EssentialGenes_List_2.txt']



all_genes_list = []
for filename in filename_list:
    file = os.path.join(path, filename)

    with open(file) as f:
        lines = f.readlines()
        print('Number of genes found in %s: %i' % (filename, (len(lines)-3)))

        for line in lines[3:]: #ASSUMES THREE HEADER LINES
            all_genes_list.append(line.strip('\n'))



del (file, filename, f, lines, line)



gene_aliases_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Yeast_Protein_Names.txt")[0]



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



unique_genes_list = list(set(all_genes_oln_list))
unique_genes_list.sort()

print('Number of unique essential genes found : %i' % len(unique_genes_list))



del (all_genes_oln_list)



#%%
if saving == True:
    save_filename = r'Cerevisiae_AllEssentialGenes_List.txt'
    
    save_file = os.path.join(path, save_filename)

    print('Creating text file with all unique genes at %s' % save_file)

    
    
    with open(save_file, 'w') as f:
        f.write('All essential genes found in lists:' + str(filename_list) + '\n')
        for gene in unique_genes_list:
            f.write(gene + '\n')

    del (gene)
