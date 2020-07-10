# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import os, sys
import timeit

dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from chromosome_and_gene_positions import gene_position
from gene_names import gene_aliases
from essential_genes_names import list_known_essentials

#%%

files_path = os.path.join(dirname, 'Data_files')
# GET POSITION GENES
gff_path = os.path.join(files_path,'Saccharomyces_cerevisiae.R64-1-1.99.gff3')
genecoordinates_dict = gene_position(gff_path) #contains all genes, essential and nonessential

# GET ALL ANNOTATED ESSENTIAL GENES
essential_path1 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_1.txt')
essential_path2 = os.path.join(files_path,'Cervisiae_EssentialGenes_List_2.txt')
essentialnames_list = list_known_essentials([essential_path1, essential_path2])

# GET ALIASES OF ALL GENES
names_path = os.path.join(files_path,'Yeast_Protein_Names.txt')
aliases_designation_dict = gene_aliases(names_path)[0]

# FOR ALL ANNOTATED ESSENTIAL GENES, DETERMINE THEIR ALIASES AND GET THE POSITION FROM GENE_POS_DICT
essentialcoordinates_dict = {}

#start_timer = timeit.default_timer()
#for gene in genecoordinates_dict:
#    if gene in essentialnames_list:
#        essentialcoordinates_dict[gene] = genecoordinates_dict.get(gene)
#    else:
#        gene_aliases_list = []
#        for key, val in aliases_designation_dict.items():
#            if gene == key or gene in val: #if gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
#                gene_aliases_list.append(key)
#                for aliases in aliases_designation_dict.get(key):
#                    gene_aliases_list.append(aliases)
#    
#        for gene_alias in gene_aliases_list:
#            if gene_alias in essentialnames_list:
#                essentialcoordinates_dict[gene_alias] = genecoordinates_dict.get(gene)
#                break
#stop_timer = timeit.default_timer()
#print('Time: %.2f seconds' %(stop_timer - start_timer))

#%%
start_timer = timeit.default_timer()
for gene in essentialnames_list:
    if gene in genecoordinates_dict:
#        print('Gene ', gene, 'in genecoordinates_dict with value ', genecoordinates_dict.get(gene))
        essentialcoordinates_dict[gene] = genecoordinates_dict.get(gene)
    else:
        gene_aliases_list = []
        for key, val in aliases_designation_dict.items():
            if gene == key or gene in val: #if gene occurs as key or in the values list in aliases_designation_dict, put all its aliases in a single list.
                gene_aliases_list.append(key)
                for aliases in aliases_designation_dict.get(key):
                    gene_aliases_list.append(aliases)

#        print('For gene ', gene, 'the following aliases are found: ', gene_aliases_list)
        for gene_alias in gene_aliases_list:
            if gene_alias in genecoordinates_dict:
#                print('The alias ', gene_alias, 'is found with value ', genecoordinates_dict.get(gene_alias))
                essentialcoordinates_dict[gene_alias] = genecoordinates_dict.get(gene_alias)
                break
#    print('')
stop_timer = timeit.default_timer()
print('Time: %.2f seconds' %(stop_timer - start_timer))

saving_list = []
for keys in essentialcoordinates_dict:
    saving_list.append(keys)