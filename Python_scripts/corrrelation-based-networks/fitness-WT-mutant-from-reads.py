 # To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'



# %% Importing libraries
from IPython import get_ipython
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.optimize as opt
from sklearn import preprocessing
from collections import defaultdict
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt


# %% Importing information per gene 
gene_info=pd.read_excel(r'chromosomal-info-per-gene.xlsx')

# %% Importing all annoatetd interactions from SGD
data_all_interactors=pd.read_excel(r'C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\Calculations\data_sgd\interaction-filtered-data.xlsx',header=0)

# %% Reading the  reads and insertions from WT and the mutant 
import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_data_wt="WT_trimmed-sorted-bam.txt"
rel_path_data_mutant="dDpl1Kan-sorted-bam.txt"


abs_path_data_wt = os.path.join(script_dir, rel_path_data_wt)
abs_path_data_mutant = os.path.join(script_dir, rel_path_data_mutant)
data_wt = pd.read_csv(abs_path_data_wt, sep="\t", header=0)
data_mutant = pd.read_csv(abs_path_data_mutant, sep="\t", header=0)


# %% Functions 
def reads2fitness (reads_data,type):
    data_reads_normalized=reads_data.copy()
    # Normalization to HO locus:
    if type == 'HO':
        data_reads_normalized['fitness-'+ type]=reads_data['Number of reads per gene']/(reads_data[reads_data['Gene name']=='HO']['Number of reads per gene'].tolist()[0])
    
 
    return data_reads_normalized

def fitness2score(gene_ref,gene_to_know, feature,data_wt,data_mutant_ref):

    info_ref=data_wt[data_wt['Gene name']==gene_ref[0]][feature].tolist()[0]


    interactions=defaultdict(dict)
    for i in np.arange(0,len(gene_to_know)):
        if len(data_wt[data_wt['Gene name']==gene_to_know[i]][feature].tolist())==0:
            info_genex=None
        else:
            info_genex=data_wt[data_wt['Gene name']==gene_to_know[i]][feature].tolist()[0]
        if len(data_mutant_ref[data_mutant_ref['Gene name']==gene_to_know[i]][feature].tolist())==0:
             info_genex_doublemutant=None
        else:
            info_genex_doublemutant=data_mutant_ref[data_mutant_ref['Gene name']==gene_to_know[i]][feature].tolist()[0]

        interactions[i]['background']=gene_ref
        interactions[i]['array-name']=gene_to_know[i]
        interactions[i]['fitness-query']=info_ref
        interactions[i]['fitness-array']=info_genex
        interactions[i]['fitness-doublemutant']=info_genex_doublemutant
        if info_genex_doublemutant==None or info_ref==None or info_genex==None:
            interactions[i]['score']=None
        else:
            interactions[i]['score']=info_genex_doublemutant-info_ref*info_genex
  
    
    return interactions

# %% Fitness from reads of WT and mutant
data_WT_normalized=reads2fitness(data_wt,type='HO')
data_mutant_normalized=reads2fitness(data_mutant,type='HO')


# %% Plotting the distribution of the fitness data
fig, axes=plt.subplots(1,2)
plt.subplots_adjust(wspace=1,right=0.5)

# Get all the data between the 1st and 3rd quartile
data=data_WT_normalized['fitness-HO']
data_iqr = data[ (data >  np.percentile(data, 25)) & (data <  np.percentile(data, 75)) ]


sns.boxplot(y='fitness-HO',data=data_WT_normalized,ax=axes[0],fliersize=0.5,color='white')
sns.stripplot(y="fitness-HO", color='black',size=2, alpha=0.1, data=data_WT_normalized,ax=axes[0])
sns.stripplot(y=data_iqr, color='green', size=2, alpha=0.2,ax=axes[0])
# axes[0].text(0.2, 0.8, str(np.round(len(data_iqr)/len(data),decimals=1)*100)+'% of the data' , horizontalalignment='center')
axes[0].set_title('WT')
# Get all the data between the 1st and 3rd quartile

data=data_mutant_normalized['fitness-HO']
data_iqr = data[ (data >  np.percentile(data, 25)) & (data <  np.percentile(data, 75)) ]

axes[0].set_ylim([0,2])
sns.boxplot(y='fitness-HO',data=data_mutant_normalized,ax=axes[1],fliersize=0.5,color='white')
sns.stripplot(y="fitness-HO", color='black',size=2, alpha=0.1, data=data_mutant_normalized,ax=axes[1])
sns.stripplot(y=data_iqr, color='green', size=2, alpha=0.2,ax=axes[1])
# axes[1].text(0.2, 0.8, str(np.round(len(data_iqr)/len(data),decimals=1)*100)+'% of the data' , horizontalalignment='center')
axes[1].set_ylim([0,2])
axes[1].set_title('mutant')


# %% evaluating the fitness score in the dpl1 mutant
interactions_dpl1=fitness2score(gene_ref=['DPL1'],gene_to_know=gene_info['gene-standard-name'],feature='fitness-HO',data_wt=data_WT_normalized,data_mutant_ref=data_mutant_normalized)


# %% Converting the dict into a dataframe 
interactions=pd.DataFrame(interactions_dpl1)
interactions_pd=interactions.T
interactions_pd.reset_index()
interactions_pd.dropna(inplace=True)
interactions_pd=interactions_pd.set_index(np.arange(0,len(interactions_pd)))


# %% Fitness plot 
fig, axes=plt.subplots(1,1)
plt.scatter(x=interactions_pd['fitness-array'],y=interactions_pd['fitness-doublemutant'],alpha=0.1)
plt.xlim([-0.05,1])
plt.xlabel('Single mutant fitness-b')
plt.ylabel('double mutant fitness-ab')
plt.ylim([-0.05,1])
x = np.linspace(0, 1)
plt.plot(x, 0.19*x,linestyle='solid',color='black');
x_masking = np.linspace(0, 0.19)
plt.plot(x, x,linestyle='solid',color='black');
plt.hlines(y=0.19,xmin=0.19,xmax=1)
plt.fill_between(x=x_masking,y1=0.19*x_masking,y2=x_masking,color='blue',alpha=0.2)
plt.fill_between(x=np.linspace(0.19,1),y2=0.19,y1=0.19*np.linspace(0.19,1),color='blue',alpha=0.2)
plt.title('Fitness related to HO locus, dpl1 fitness=0.19')
#plt.savefig('dpl1_data_using_reads_normalized_with_HO_as_fitness.png',format='png',dpi=300,transparent=True)




# %% Selecting columns to measure the correlations 
values2corr=interactions_pd.loc[:,['fitness-array','fitness-doublemutant','score']]
values2corr=values2corr.astype(float)
values2corr.reset_index()
values2corr=values2corr.set_index(np.arange(0,len(values2corr)))



# %% Row- wise correlation with data from dpl1 

corr_with_a_row=values2corr.corrwith(values2corr.iloc[108,:],axis=1,method='pearson')
corr_with_a_row_pd=pd.DataFrame(corr_with_a_row,columns=['corr-values'])
corr_with_a_row_pd.sort_values(by='corr-values',inplace=True)

interactions_pd.set_index=np.arange(0,len(interactions_pd))
index_corr=corr_with_a_row_pd.index
corr_with_a_row_pd['genes']=interactions_pd.ix[index_corr,'array-name']
corr_with_a_row_pd.fillna(0,inplace=True)


dpl1_interactors_data=data_all_interactors[data_all_interactors['Standard_Gene_Name_(Bait)']=='DPL1']


# %% Inspecting if the estimated correlations make sense with the existing info for Dpln1

for i in dpl1_interactors_data.iloc[:,1].unique():
    tmp=corr_with_a_row_pd[corr_with_a_row_pd['genes']==i]
    

    if len(tmp)!=0:

        tmp_2=data_all_interactors[data_all_interactors.iloc[:,1]==i]
        type_int=tmp_2[data_all_interactors.columns[2]]
        

        corr_with_a_row_pd.loc[tmp.index,'interactor']=True

        if type_int.value_counts()[:1].index.tolist()== ['Negative Genetic']:
            corr_with_a_row_pd.loc[tmp.index,'type']=-1
        elif type_int.value_counts()[:1].index.tolist()== ['Synthetic Lethality']:
            corr_with_a_row_pd.loc[tmp.index,'type']=-1
        elif type_int.value_counts()[:1].index.tolist()== ['Positive Genetic']:
            corr_with_a_row_pd.loc[tmp.index,'type']=1
        elif type_int.value_counts()[:1].index.tolist()== ['Synthetic Rescue']:
            corr_with_a_row_pd.loc[tmp.index,'type']=1
        else:
            corr_with_a_row_pd.loc[tmp.index,'type']=0.5
        corr_with_a_row_pd.loc[tmp.index,'frequency-of-the-true-type']= type_int.value_counts().max()/np.sum(type_int.value_counts())


    else:
        corr_with_a_row_pd.loc[tmp.index,'interactor']=False
        corr_with_a_row_pd.loc[tmp.index,'type']=False


corr_with_a_row_pd.fillna(0,inplace=True)  



# %% Heatmap of correlations


plt.figure(figsize=(5,5))
sns.heatmap(data=corr_with_a_row_pd[corr_with_a_row_pd['interactor']==True].corr(),center=0,vmin=-1,cmap='PiYG')



# %% Re-structuring the correlation values to discrete values
df=corr_with_a_row_pd.copy()
corr_with_a_row_pd['corr-values-group']=pd.cut(df['corr-values'], bins=5,labels=np.arange(-1,1.5,0.5), right=False)


# %% [markdown]
# ### See this https://towardsdatascience.com/visualising-stocks-correlations-with-networkx-88f2ee25362e

# %% Categoryzing new genes and existing genes that the algorithm give

for i in np.arange(0,len(corr_with_a_row_pd)):
    if corr_with_a_row_pd.loc[i,'type']== corr_with_a_row_pd.loc[i,'corr-values-group']:
        corr_with_a_row_pd.loc[i,'genes-classified']=corr_with_a_row_pd.loc[i,'genes']
    else:
        corr_with_a_row_pd.loc[i,'genes-classified']=corr_with_a_row_pd.loc[i,'genes'] + '-new'


# %% Plot to see how close is the type and the correlation values

plt.subplots(1,1)
data_dpl1=corr_with_a_row_pd[corr_with_a_row_pd['interactor']==1]
sns.scatterplot(x='corr-values-group',y='type',data=data_dpl1)


# %% Drawing the network 
# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
 

 
# Build your graph

G_estimated=nx.from_pandas_edgelist(df=corr_with_a_row_pd[corr_with_a_row_pd['interactor']==1],source='corr-values-group',target='genes-classified')

# degree of the plot
d_estimated = nx.degree(G_estimated)

#####creates list of nodes and a list their degrees that will be used later for their sizes
nodelist_estimated, node_sizes_estimated = zip(*d_estimated)

color_map = []
for node in G_estimated:
    if node == 1:
        color_map.append('green')
    elif node==-1:
        color_map.append('red')
    else: 
        color_map.append('gray')    
        
plt.figure(figsize=(5,5))
nx.draw(G_estimated,pos=nx.spring_layout(G_estimated), with_labels=True,node_color=color_map, node_size=tuple([x**3 for x in node_sizes_estimated]), edge_color='black', linewidths=1, font_size=10);

plt.savefig('corr-based-net-dpl1-test.png',dpi=300,format='png',transparent=False)






