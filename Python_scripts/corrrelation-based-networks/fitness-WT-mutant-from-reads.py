 # To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'



# %% Importing libraries
from IPython import get_ipython
import pandas as pd
import numpy as np
import seaborn as sns
from collections import defaultdict
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt


# %% Importing information per gene 
gene_info=pd.read_excel(r'C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\SATAY\src(source-code)\LaanLab-SATAY-DataAnalysis\Python_scripts\corrrelation-based-networks\chromosomal-info-per-gene.xlsx')

# %% Importing all annotated interactions from SGD
data_all_interactors=pd.read_excel(r'C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\Calculations\data_sgd\interaction-filtered-data.xlsx',header=0)

# %% Reading the  reads and insertions from WT and the mutant 
import os

script_dir = os.path.dirname(os.path.abspath('__file__')) #<-- absolute dir the script is in
rel_path_data_wt="WT-normalize-reads-all-chromosome.xlsx"
rel_path_data_mutant="dpl1-normalize-reads-all-chromosome.xlsx"
#%% Reading the datafiles 

abs_path_data_wt = os.path.join(script_dir,'Python_scripts','fitness_from SATAY_codes','output_data_files',rel_path_data_wt)
abs_path_data_mutant = os.path.join(script_dir,'Python_scripts','fitness_from SATAY_codes','output_data_files',rel_path_data_mutant)
# data_wt = pd.read_csv(abs_path_data_wt, sep="\t", header=0)
# data_mutant = pd.read_csv(abs_path_data_mutant, sep="\t", header=0)
data_wt=pd.read_excel(abs_path_data_wt,index_col='Unnamed: 1')
data_wt.fillna('',inplace=True)
data_wt.rename(columns={'Unnamed: 0':'chrom'},inplace=True)

data_mutant=pd.read_excel(abs_path_data_mutant,index_col='Unnamed: 1')
data_mutant.fillna('',inplace=True)
data_mutant.rename(columns={'Unnamed: 0':'chrom'},inplace=True)



# %% Functions 
def reads2fitness (reads_data,type,feature):
    """
    

    Parameters
    ----------
    reads_data : dataframe 
        where the data is 
    type : string
        Name of the gene to normalize with , e.g . 'HO'
    feature : string
        Name of the column to normalize with
        examples: 'Nreadsperbp_central80p_normalized', 'Nreadsperbp_normalized'

    Returns
    -------
    data_reads_normalized : dataframe
        same dataframe as the input with an extra column with the fitness values 

    """
    data_reads_normalized=reads_data.copy()
    # Normalization to HO locus:
    if type == 'HO':
        data_reads_normalized['fitness-'+ type]=reads_data[feature]/(reads_data[reads_data['Standard_name']==type][feature].tolist()[0])
    
 
    return data_reads_normalized

def fitness2score(gene_ref,gene_to_know, feature,data_wt,data_mutant_ref):

    info_ref=data_wt[data_wt['Standard_name']==gene_ref[0]][feature].tolist()[0]


    interactions=defaultdict(dict)
    for i in np.arange(0,len(gene_to_know)):
        if len(data_wt[data_wt['Standard_name']==gene_to_know[i]][feature].tolist())==0:
            info_genex=None
        else:
            info_genex=data_wt[data_wt['Standard_name']==gene_to_know[i]][feature].tolist()[0]
        if len(data_mutant_ref[data_mutant_ref['Standard_name']==gene_to_know[i]][feature].tolist())==0:
             info_genex_doublemutant=None
        else:
            info_genex_doublemutant=data_mutant_ref[data_mutant_ref['Standard_name']==gene_to_know[i]][feature].tolist()[0]

        interactions[i]['background']=gene_ref
        interactions[i]['array-name']=gene_to_know[i]
        interactions[i]['fitness-query']=info_ref
        interactions[i]['fitness-array']=info_genex
        interactions[i]['fitness-doublemutant']=info_genex_doublemutant
        if info_genex_doublemutant==None or info_ref==None or info_genex==None:
            interactions[i]['score']=None
        else:
            interactions[i]['score']=info_genex_doublemutant-(info_ref*info_genex)
  
    
    return interactions

# %% Fitness from reads of WT and mutant
data_WT_normalized=reads2fitness(data_wt,type='HO',feature='Nreads')
data_mutant_normalized=reads2fitness(data_mutant,type='HO',feature='Nreads')


# %% Plotting the distribution of the fitness data
fig, axes=plt.subplots(1,2)
plt.subplots_adjust(wspace=1,right=0.8)

# Get all the data between the 1st and 3rd quartile
data=data_WT_normalized['fitness-HO']
data_iqr = data[ (data >  np.percentile(data, 25)) & (data <  np.percentile(data, 75)) ]


sns.boxplot(y='fitness-HO',data=data_WT_normalized,ax=axes[0],fliersize=0.5,color='white')
sns.stripplot(y="fitness-HO", color='black',size=2, alpha=0.1, data=data_WT_normalized,ax=axes[0])
#sns.stripplot(y=data_iqr, color='green', size=2, alpha=0.1,ax=axes[0])
# axes[0].text(0.2, 0.8, str(np.round(len(data_iqr)/len(data),decimals=1)*100)+'% of the data' , horizontalalignment='center')
axes[0].set_title('WT')
# Get all the data between the 1st and 3rd quartile

data=data_mutant_normalized['fitness-HO']
data_iqr = data[ (data >  np.percentile(data, 25)) & (data <  np.percentile(data, 75)) ]

axes[0].set_ylim([0,20])
sns.boxplot(y='fitness-HO',data=data_mutant_normalized,ax=axes[1],fliersize=0.5,color='white',saturation=0.5)
sns.stripplot(y="fitness-HO", color='black',size=2, alpha=0.1, data=data_mutant_normalized,ax=axes[1])
#sns.stripplot(y=data_iqr, color='green', size=2, alpha=0.1,ax=axes[1])
# axes[1].text(0.2, 0.8, str(np.round(len(data_iqr)/len(data),decimals=1)*100)+'% of the data' , horizontalalignment='center')
axes[1].set_ylim([0,20])
axes[1].set_title('mutant-dpl1d')
#%%
fig1, axes=plt.subplots(1,1)
plt.subplots_adjust(wspace=1,right=0.8)
plt.scatter(x=data_WT_normalized['fitness-HO'],y=data_mutant_normalized['fitness-HO'],color='black',alpha=0.2)
plt.xlim([0,50])
plt.ylim([0,50])
plt.xlabel('fitness-HO-WT')
plt.ylabel('fitness-HO-dpl1d')
plt.plot(np.linspace(0,100),np.linspace(0,100))

#%% saving the figure

fig.savefig('fitness-based-on-HO-locus.png',dpi=300,transparency=True,format='png')
fig1.savefig('scatter-plot-fitness-HO-dpl1-vs-WT.png',dpi=300,transparency=True,format='png')

# %% evaluating the fitness score in the dpl1 mutant
interactions_dpl1=fitness2score(gene_ref=['DPL1'],gene_to_know=gene_info['gene-standard-name'],feature='fitness-HO',data_wt=data_WT_normalized,data_mutant_ref=data_mutant_normalized)


# %% Converting the dict into a dataframe 
interactions=pd.DataFrame(interactions_dpl1)
interactions_pd=interactions.T
interactions_pd.reset_index()
interactions_pd.dropna(inplace=True)
interactions_pd=interactions_pd.set_index(np.arange(0,len(interactions_pd)))

#%%
dpl1_fitness=interactions_pd[interactions_pd['array-name']=='DPL1']['fitness-query'].tolist()[0]
# %% Fitness plot 
max=10
min=-1
fig, axes=plt.subplots(1,1)

plt.scatter(x=interactions_pd['fitness-array'],y=interactions_pd['fitness-doublemutant'],alpha=0.1)
plt.xlim([min,max])
plt.xlabel('Single mutant fitness-b')
plt.ylabel('double mutant fitness-ab')
plt.ylim([min,max])
x = np.linspace(min, max)

x_masking = np.linspace(min, dpl1_fitness)
plt.plot(x, x,linestyle='solid',color='gray');
#plt.hlines(y=dpl1_fitness,xmin=dpl1_fitness,xmax=max)
#plt.fill_between(x=x_masking,y1=dpl1_fitness*x_masking,y2=x_masking,color='blue',alpha=0.2)
#plt.fill_between(x=np.linspace(dpl1_fitness,max),y2=dpl1_fitness,y1=dpl1_fitness*np.linspace(dpl1_fitness,max),color='blue',alpha=0.2)
plt.title('Fitness related to HO locus, dpl1 fitness='+ str(np.round(dpl1_fitness,2)))

plt.plot(x, dpl1_fitness*x,linestyle='solid',color='gray');


trianglex=[0,5,5,0]
triangley=[0,dpl1_fitness*5,5,0]

#plt.fill(trianglex, triangley,alpha=0.3,color='gray')
#%% saving the figure 
fig.savefig('dpl1_data_using_reads_normalized_with_HO_as_fitness.png',format='png',dpi=300,transparent=True)

# %% Selecting columns to measure the correlations 
values2corr=interactions_pd.loc[:,['fitness-array','fitness-doublemutant','score']]
values2corr=values2corr.astype(float)
values2corr.reset_index()
values2corr=values2corr.set_index(np.arange(0,len(values2corr)))



# %% Row- wise correlation with data from dpl1 
# - try to do a loop to do different corrwith based on different genes and then average and see the std of the different correlations 

corr_with_a_row=values2corr.corrwith(values2corr.iloc[997,:],axis=1,method='pearson') # dpl1 row
#corr_with_a_row=values2corr.corrwith(values2corr.iloc[interactions_pd[interactions_pd['array-name']=='AIM22'].index,:],axis=1,method='pearson')
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

        # tmp_2=data_all_interactors[data_all_interactors.iloc[:,1]==i]
        # type_int=tmp_2[data_all_interactors.columns[2]]
        

        corr_with_a_row_pd.loc[tmp.index,'interactor']=True
        
        type_int=dpl1_interactors_data[dpl1_interactors_data.iloc[:,1]==i]['Experiment_Type_(Required) ']

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



# Heatmap of correlations


plt.figure(figsize=(5,5))
sns.heatmap(data=corr_with_a_row_pd[corr_with_a_row_pd['interactor']==True].corr(),center=0,vmin=-1,cmap='PiYG')



# %% Re-structuring the correlation values to discrete values
df=corr_with_a_row_pd.copy()
corr_with_a_row_pd['corr-values-group']=pd.cut(df['corr-values'], bins=5,labels=np.arange(-1,1.5,0.5), right=False)




#  Categoryzing new genes and existing genes that the algorithm give

for i in np.arange(0,len(corr_with_a_row_pd)):
    if corr_with_a_row_pd.loc[i,'type']== corr_with_a_row_pd.loc[i,'corr-values-group']:
        corr_with_a_row_pd.loc[i,'genes-classified']=corr_with_a_row_pd.loc[i,'genes']
    else:
        corr_with_a_row_pd.loc[i,'genes-classified']=corr_with_a_row_pd.loc[i,'genes'] + '-new'


# Plot to see how close is the type and the correlation values

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


#%% dendogram

from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

# ytdist = np.array([662., 877., 255., 412., 996., 295., 468., 268.,
#                    400., 754., 564., 138., 219., 869., 669.])

#ytdist=np.array(corr_with_a_row_pd['corr-values'])

ytdist=np.array(values2corr['score'][0:10])
Z = hierarchy.linkage(np.reshape(ytdist, (len(ytdist), 1)))
plt.figure(figsize=(8,8))
dn = hierarchy.dendrogram(Z)

# %% 1D heatmap
import matplotlib.pyplot as plt
import numpy as np; np.random.seed(1)
plt.rcParams["figure.figsize"] = 5,2

x = np.linspace(0,len(values2corr['score']))
x=np.arange(0,len(values2corr['score']))
# y = np.cumsum(np.random.randn(50))+6
y=values2corr['score']
 

fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)

#extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
ax.imshow(y[np.newaxis,:], cmap="Blues", aspect='auto',vmin=-1,vmax=2)
ax.set_yticks([])
# ax.set_xlim(extent[0], extent[1])


ax2.plot(x,y[x],'*')
ax2.set_ylim([-5,10])

plt.tight_layout()
plt.show()
