# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 10:34:28 2020

@author: Thomas
"""

from intermine.webservice import Service
import numpy as np
import scipy.io
import seaborn as sns
from scipy import stats, optimize, interpolate
import pandas as pd
from collections import defaultdict 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from scipy import stats
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os, fnmatch

import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
#rel_path_domains="datasets/proteins-domains-from-Pfam.xlsx"

# abs_file_path_domains = os.path.join(script_dir, rel_path_domains)

# os.chdir('../') #<-- for binder os.chdir('../')
# my_path_domains=abs_file_path_domains
#data_domains=pd.read_excel(my_path_domains,header=0,index_col='Unnamed: 0')
#data_domains=data_domains.dropna()

#data_domains.head()

pathways=[ 'galactose degradation','phospholipid biosynthesis', "ethanol degradation"]
biological_processes_slim=["cell budding","lipid binding","cytokinesis"]
biological_processes_goterm=["cell morphogenesis involved in conjugation with cellular fusion","budding cell apical bud growth","establishment of cell polarity","cytoskeleton organization"]
def from_go_to_genes(go,label):
    #label=["GOTerm" or "GOSlimTerm"]
    service = Service('https://yeastmine.yeastgenome.org/yeastmine/service', token = 's1Z7G6d1KdH91d28kd70')
    query = service.new_query("Gene")
    query.add_constraint("goAnnotation.ontologyTerm.parents", label)
    query.add_view(
        "symbol", "goAnnotation.evidence.code.annotType",
        "goAnnotation.ontologyTerm.parents.name"
    )
    query.add_constraint("goAnnotation.qualifier", "!=", "NOT", code = "C")
    query.add_constraint("goAnnotation.qualifier", "IS NULL", code = "D")
    query.add_constraint("goAnnotation.evidence.code.annotType", "=", "manually curated", code = "F")
    query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", go, code = "G")
    query.set_logic("(C or D) and F and G")

    data_toy=defaultdict(dict)

    for row,counter in zip(query.rows(),np.arange(0,len(query.rows()))):

        data_toy['gene-name'][counter]=row["symbol"]
        data_toy['evidence'][counter]=row["goAnnotation.evidence.code.annotType"]
        data_toy['annotation'][counter]=row["goAnnotation.ontologyTerm.parents.name"]


    data_toy_pd=pd.DataFrame(data_toy)
    data=data_toy_pd.drop_duplicates()

    data.index=np.arange(0,len(data))
    return data

#%%

go=biological_processes_goterm[2]
data=from_go_to_genes(go,label='GOTerm')
#%%
big_data=[]


for i in np.arange(0,len(biological_processes_goterm)):
    go=biological_processes_goterm[i]
    data=from_go_to_genes(go,label='GOTerm')
    big_data.append(data['gene-name'].tolist())
    
flat_list = [item for sublist in big_data for item in sublist]
res = [i for i in flat_list if i] 
flat_list_unique=np.unique(res)

# Selecting the meaningful columns in the respective dataset
domain_id_list=data_domains['domain-name']
#query_gene=data['gene-name']
query_gene=flat_list_unique

# Initialising the arrays
protein_a_list=[]

population = np.arange(0,len(query_gene))


for m in population:

    protein_a=data_domains[data_domains['name']==query_gene[m]]
    protein_a_list.append(protein_a['domain-name'].tolist())
    


def remove_empty_domains(protein_list_search):
    index=[]
    for i in np.arange(0,len(protein_list_search)):
        if protein_list_search[i]==[]:
            index.append(i) ## index of empty values for the protein_a_list meaning they dont have any annotated domain

    y=[x for x in np.arange(0,len(protein_list_search)) if x not in index] # a list with non empty values from protein_a list

    protein_list_search_new=[]
    
    for i in y:
        protein_list_search_new.append(protein_list_search[i])
        
    return protein_list_search_new

## evaluating the function

protein_a_list_new=remove_empty_domains(protein_a_list)

print('The empty domain in the data were:', len(protein_a_list)-len(protein_a_list_new), 'out of', len(protein_a_list),'domains')

protein_module=[]

for i in np.arange(0,len(protein_a_list_new)):
    protein_module.append(np.unique(protein_a_list_new[i]))# just taking the unique domains per protein in the module 
protein_module=np.concatenate(protein_module, axis=0 )
unique_domains,index_domains,index_counts=np.unique(protein_module,return_index=True,return_counts=True) # counting how many proteins have the same domain per module 

stats_module=pd.DataFrame()
stats_module['unique-domains']=unique_domains
stats_module['counts-each-domains']=index_counts
stats_module['index-for-original-domains-array']=index_domains
stats_module['ratio-of-sharing']=(index_counts-1)/len(unique_domains)

stats_module.set_index=np.arange(0,len(stats_module))