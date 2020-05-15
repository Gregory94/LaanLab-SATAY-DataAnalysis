'''
This script reads:
- 'Data-fitness.xlsx'
- 'transposon-density-per-gene-benoit.xlsx'.
It also requires:
- r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_1.txt"
- r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_2.txt"
- r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt"
and requires the python modules
- essential_genes_names.py
- gene_names.py

This script creates a scatterplot for transposon_density vs fitness-SGA.
All genes that are used to create the plot are searched for is they are annotated as essential genes.
This is also done for any potential aliases the genes might have.

This code is a part from 'fitness_vs_transposon_density.py' created by Leila.
'''

#import numpy as np
#import scipy.io
import seaborn as sns
#from scipy import stats, optimize, interpolate
import pandas as pd
from collections import defaultdict 
#import math
import matplotlib.pyplot as plt
#from scipy.stats import norm, lognorm
#from scipy import stats
#import matplotlib.cm as cm
#import matplotlib.mlab as mlab
#import os, fnmatch
import sys
#SET PATH TO PYTHON MODULES (r"N:\tnw\BN\LL\Shared\Gregory\Python\Python_modules")
sys.path.insert(1,r"C:\Users\gregoryvanbeek\Documents\Python Scripts\Python_modules")
from essential_genes_names import list_known_essentials
from gene_names import gene_aliases

## published data in the cellmap.org 2016 
data_fitness_sga=pd.read_excel(r"C:\Users\gregoryvanbeek\Documents\Gene_Database\Fitness_Data\Data-fitness.xlsx",header=0,sheet_name='NxN')

data_fitness_sga.columns=['query-strain','query-allele-name','array-strain','array-allele-name','array-type','score','p-value','query-fitness','array-fitness','double-fitness','double-fitness-std']

## data saved from the df dataframe of statistics_per_gene.ipynb
data_transposon=pd.read_excel(r"C:\Users\gregoryvanbeek\Documents\Gene_Database\transposon-density-per-gene-benoit.xlsx",header=0)
data_transposon=data_transposon.drop(['Unnamed: 0'],axis=1)
data_transposon=data_transposon.apply(lambda x: x.astype(str).str.lower()) # make everything lowercase

#transform back these columns into float variables
data_transposon['Transposon_density_per_gene']=data_transposon['Transposon_density_per_gene'].apply(lambda x: (float(x)))
data_transposon['Read_density_per_gene']=data_transposon['Read_density_per_gene'].apply(lambda x: (float(x)))




df_transposon_fitness=defaultdict(dict)

for i in data_transposon['Gene_name']:
    if len(data_fitness_sga[data_fitness_sga['query-allele-name']==i])!=0:
        df_transposon_fitness['gene-merge'][i]=i
        df_transposon_fitness['transposon-density'][i]=float(data_transposon[data_transposon['Gene_name']==i]['Transposon_density_per_gene'].tolist()[0])
        df_transposon_fitness['fitness-SGA'][i]=data_fitness_sga[data_fitness_sga['query-allele-name']==i]['query-fitness'].tolist()[0]
        if data_transposon[data_transposon['Gene_name']==i]['Essential_gene'].tolist()==['true']:
            df_transposon_fitness['essential-by-benoit']='true'
        else :
            df_transposon_fitness['essential-by-benoit']='false'
    else :
        if len(data_fitness_sga[data_fitness_sga['array-allele-name']==i])!=0:
            df_transposon_fitness['gene-merge'][i]=i
            df_transposon_fitness['transposon-density'][i]=float(data_transposon[data_transposon['Gene_name']==i]['Transposon_density_per_gene'].tolist()[0])
            df_transposon_fitness['fitness-SGA'][i]=data_fitness_sga[data_fitness_sga['array-allele-name']==i]['array-fitness'].tolist()[0]
            if data_transposon[data_transposon['Gene_name']==i]['Essential_gene'].tolist()==['true']:
                df_transposon_fitness['essential-by-benoit']='true'
            else :
                df_transposon_fitness['essential-by-benoit']='false'

df_transposon_fitness=pd.DataFrame(df_transposon_fitness)
df_transposon_fitness_filtered=df_transposon_fitness[df_transposon_fitness['fitness-SGA']<1]



#USE FUNCTIONS TO GET NAMES OF ESSENTIAL GENES. FILES USED ARE:
#r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_1.txt"
#r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_2.txt"
#r"N:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt"
known_essential_gene_list = list_known_essentials([r"C:\Users\gregoryvanbeek\Documents\Gene_Database\Cervisiae_EssentialGenes_List_1.txt",r"C:\Users\gregoryvanbeek\Documents\Gene_Database\Cervisiae_EssentialGenes_List_2.txt"])
aliases_designation_dict, aliases_sgd_dict, aliases_swissprot_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\Gene_Database\Yeast_Protein_Names.txt")

essential_counter = 0
for n in range(df_transposon_fitness_filtered.shape[0]):
    gene_inDF = df_transposon_fitness_filtered['gene-merge'][n].upper()
    if gene_inDF in known_essential_gene_list:
        df_transposon_fitness_filtered['essential-by-benoit'][n] = 'true'
        print(gene_inDF, ' is annotated as essential')
        essential_counter += 1
    elif gene_inDF in aliases_designation_dict:
        for alias in aliases_designation_dict.get(gene_inDF):
            if alias in known_essential_gene_list:
                df_transposon_fitness_filtered['essential-by-benoit'][n] = 'true'
                print(gene_inDF, ' with alias ',alias, ' is annotated as essential')
                essential_counter += 1
    else:
        for key, val in aliases_designation_dict.items():
            if gene_inDF in val:
                for alias in val:
                    if alias in known_essential_gene_list:
                        df_transposon_fitness_filtered['essential-by-benoit'][n] = 'true'
                        print(gene_inDF, ' with alias ',alias, ' is annotated as essential')
                        essential_counter += 1
print('Number of essential genes found = ',essential_counter)


df_transposon_fitness_filtered_renamed = df_transposon_fitness_filtered.rename(columns={'essential-by-benoit': 'annotated-essentials'})
df_transposon_fitness_filtered_ordered = df_transposon_fitness_filtered_renamed.sort_values(by=['annotated-essentials'], ascending=True)
sns.scatterplot(x='transposon-density' ,y='fitness-SGA',hue='annotated-essentials', data=df_transposon_fitness_filtered_ordered,palette=['green','red'],alpha=0.6, linewidth=0)
#plt.savefig('transposon-density-vs-fitness-NxN-scatter.png',dpi=300,format='png',transparent=False)
#sns.jointplot("transposon-density", "fitness-SGA", data=df_transposon_fitness_filtered, height=6, ratio=3, color='essential-by-benoit', xlim=[0,0.5],ylim=[0,1.0],kind='scatter',alpha=0.3)



