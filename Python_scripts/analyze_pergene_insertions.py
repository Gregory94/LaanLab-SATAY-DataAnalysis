# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 08:25:34 2020

@author: gregoryvanbeek

THIS SCRIPT ANALYSE THE INSERTION LOCATIONS PER GENE AS STORED IN _PERGENE_INSERTIONS.TXT AND PERESSENTIAL_INSERTIONS.TXT
"""

#%%
import os, sys
import copy
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.cbook import boxplot_stats


dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(dirname,'python_modules'))
from gene_names import gene_aliases

#%%
def tninserts_analysis():
    '''
    IMPORTANT VARIABLES IN THIS FUNCTION:

    gene_position_dict: position of all genes
    gene_inserts_dict: insertion locations of transposons for all genes
    gene_reads_dict: number of reads for all tn inserts for all genes

    essential_position_dict: only for essential genes.
    essential_inserts_dict: only for essential genes.
    essential_reads_dict: only for essential genes.

    nonessential_position_dict: only for nonessential genes.
    nonessential_inserts_dict: only for nonessential genes.
    nonessential_reads_dict: only for nonessential genes.
    
    df: dataframe to store all information for analysis
    '''
#%% READ FILE AND PUT ALL VALUES IN DICTIONARIES. DO NOT CHANGE THIS SECTION.
    filepath = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out"
    filename = "ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"
#    filepath = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt2_testfolder\WT2_dataset_analysis_temp202008051429_new2"
#    filename = r"E-MTAB-4885.WT2.bam_pergene_insertions.txt"
    datafile = os.path.join(filepath, filename)

    with open(datafile) as f:
        lines = f.readlines()


    gene_position_dict = {}
    gene_inserts_dict = {}
    gene_reads_dict = {}
    
    gene_inserts_distance_dict = {} #distance between subsequent inserts
    gene_inserts_trunc_dict = {} #inserts in the gene where 10% of the edges is truncated (so, only the center part of the gene is considered).
    gene_reads_trunc_dict = {} #reads in the gene where 10% of the edges is truncated (so, only the center part of the gene is considered).
    for line in lines[1:]:
        line_split = line.strip('\n').split('\t')
        genename = line_split[0]
        gene_chrom = line_split[1]
        gene_start = int(line_split[2])
        gene_end = int(line_split[3])


        gene_position_dict[genename] = [gene_chrom, gene_start, gene_end]


        geneinserts_str = line_split[4].strip('[]')
        if not geneinserts_str == '':
            geneinserts_list = [int(ins) for ins in geneinserts_str.split(',')]
        else:
            geneinserts_list = []
        gene_inserts_dict[genename] = geneinserts_list


        ins_list = []
        ins_indx_list = []
        for ins in geneinserts_list: #GET INSERTIONS THAT ARE MORE THAN 10% OF LENGTH GENE AWAY FROM THE EDGES OF THE GENE.
            l = gene_end - gene_start
            if (gene_start + 0.1*l) < ins < (gene_end - 0.1*l):
                ins_list.append(ins)
                ins_indx_list.append(geneinserts_list.index(ins))
        gene_inserts_trunc_dict[genename] = ins_list


        if not len(geneinserts_list) < 2:
            d = []
            for i in range(1,len(geneinserts_list)): #DISTANCES BETWEEN SUBSEQUENT INSERTS
                d.append(geneinserts_list[i] - geneinserts_list[i-1])
            gene_inserts_distance_dict[genename] = d
        elif len(geneinserts_list) == 1:
            gene_inserts_distance_dict[genename] = [0]#[0] #only one insert
        else:
            gene_inserts_distance_dict[genename] = [0]#[-1] #no insert


        genereads_str = line_split[5].strip('[]')
        if not genereads_str == '':
            genereads_list = [int(read) for read in genereads_str.split(',')]
        else:
            genereads_list = []
        gene_reads_dict[genename] = genereads_list
        gene_reads_trunc_dict[genename] = [genereads_list[i] for i in ins_indx_list]

        if len(geneinserts_list) != len(genereads_list):
            print('WARNING: %s has different number of reads compared with the number of inserts' % genename )



    del (datafile, lines, line, line_split, genename, gene_chrom, gene_start, gene_end, geneinserts_str, geneinserts_list, genereads_str, genereads_list, i, d, ins, ins_list, ins_indx_list, l)
    #remains: gene_inserts_dict, gene_position_dict, gene_reads_dict


#%% GET ANNOTATED ESSENTIAL GENES. DO NOT CHANGE THIS SECTION.
    essentialsfile = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_AllEssentialGenes_List.txt"

    with open(essentialsfile) as f:
        lines = f.readlines()


    aliases_dict = gene_aliases(r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Yeast_Protein_Names.txt")[0]


    essential_position_dict = {}
    essential_inserts_dict = {}
    essential_reads_dict = {}

    nonessential_position_dict = copy.deepcopy(gene_position_dict)
    nonessential_inserts_dict = copy.deepcopy(gene_inserts_dict)
    nonessential_reads_dict = copy.deepcopy(gene_reads_dict)


    for line in lines[1:]:
        essential = line.strip('\n')


        essentiality = 'nonessential'


        if essential in gene_position_dict:
            essentiality = 'essential'
            alias = essential
        else:
            for alias in aliases_dict.get(essential):
                if alias in gene_position_dict:
                    essentiality = 'essential'
                    break


        if essentiality == 'essential':
            essential_position_dict[alias] = gene_position_dict.get(alias)
            essential_inserts_dict[alias] = gene_inserts_dict.get(alias)
            essential_reads_dict[alias] = gene_reads_dict.get(alias)
            
            del nonessential_position_dict[alias]
            del nonessential_inserts_dict[alias]
            del nonessential_reads_dict[alias]


    del (essentialsfile, lines, line, aliases_dict, essential, essentiality, alias)
    #remain: essential_position_dict, essential_inserts_dict, essential_reads_dict, nonessential_position_dict, nonessential_inserts_dict, nonessential_reads_dict



#%% CREATE DATAFRAME FOR ALL GENES. ADD STATISTICS HERE
    genename_list = []
    essentiality_list = []
    N_inserts_notnormalized_list = []
    N_reads_notnormalized_list = []
    N_inserts_list = []
    N_inserts_trunc_list = []
    N_reads_list = []
    N_reads_trunc_list = []
    N_reads_per_tn_full_list = []
    N_reads_per_tn_trunc_list = []
    distance_max_inserts_list = []
    gene_id_list = []
    i = 0
    for gene in gene_position_dict:
        genename_list.append(gene) #GENENAME LIST
        gene_id_list.append(i)
        i += 1


        if gene in essential_position_dict: #ESSENTIALITY_LIST
            essentiality_list.append(True)
        elif gene in nonessential_position_dict:
            essentiality_list.append(False)
        else:
            print('WARNING: %s not found.' % gene)


        l_gene = gene_position_dict.get(gene)[2] - gene_position_dict.get(gene)[1]


        N_inserts_notnormalized_list.append(len(gene_inserts_dict.get(gene)))
        N_reads_notnormalized_list.append(sum(gene_reads_dict.get(gene)))


        N_inserts_list.append(len(gene_inserts_dict.get(gene)) / l_gene) #N_INSERTS_LIST (NUMBER OF INSERTIONS)
        N_inserts_trunc_list.append(len(gene_inserts_trunc_dict.get(gene)) / l_gene) #N_INSERTS_CENTER_LIST (NUMBER OF INSERTIONS IN THE GENE WHERE 10% OF THE GENE LENGTH IS TRUNCATED)

        N_reads_list.append(sum(gene_reads_dict.get(gene)) / l_gene) #N_READS_LIST (TOTAL NUMBER OF READS IN GENE)
        N_reads_trunc_list.append(sum(gene_reads_trunc_dict.get(gene)) / l_gene)

        if not len(gene_inserts_dict.get(gene)) == 0:
            N_reads_per_tn_full_list.append(sum(gene_reads_dict.get(gene)) / len(gene_inserts_dict.get(gene)))
        else:
            N_reads_per_tn_full_list.append(0)

        if not gene_inserts_trunc_dict.get(gene) == []:
            N_reads_per_tn_trunc_list.append(sum(gene_reads_trunc_dict.get(gene)) / len(gene_inserts_trunc_dict.get(gene)))
        else:
            N_reads_per_tn_trunc_list.append(0)

        distance_max_inserts_list.append(np.nanmax(gene_inserts_distance_dict.get(gene)) / l_gene) #DISTANCE_MAX_INSERTS_LIST (LARGEST DISTANCE BETWEEN SUBSEQUENT INSERTIONS NORMALIZED TO GENE LENGTH)



    allgenes = {'Gene_Name': genename_list,
                'Gene_ID': gene_id_list,
                'Essentiality': essentiality_list,
                'Insertions_NotNorm_Full_Gene': N_inserts_notnormalized_list,
                'Reads_NotNorm_Full_Gene': N_reads_notnormalized_list,
                'Number_Insertions_Full_Gene': N_inserts_list,
                'Number_Insertions_Truncated_Gene': N_inserts_trunc_list,
                'Number_Reads_Full_Gene': N_reads_list,
                'Number_Reads_Truncated_Gene': N_reads_trunc_list,
                'Reads_per_Transposon_full_Gene': N_reads_per_tn_full_list,
                'Reads_per_Transposon_Truncated_Gene': N_reads_per_tn_trunc_list,
                'Max_Insertion_Distance': distance_max_inserts_list}


    df = pd.DataFrame(allgenes, columns = [column_name for column_name in allgenes])


    del (gene, genename_list, gene_id_list, i, l_gene, essentiality_list, N_inserts_list, N_inserts_trunc_list, N_reads_per_tn_trunc_list, distance_max_inserts_list, N_reads_per_tn_full_list, N_reads_list, N_reads_trunc_list, allgenes)


#%%TEST GRAPH
    sns.set(style="whitegrid")

    #POTENTIALLY USEFUL; NUMBER OF INSERTIONS IN THE ENTIRE GENE.
    ax1 = sns.boxplot(x='Essentiality',y='Number_Insertions_Full_Gene',data=df)
    ax1.set_ylim(-0.001,0.04)
    ax1.set_xlabel('Annotated essential', fontsize=14)
    ax1.set_ylabel('Number of insertions (normalized to gene length)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)


    #USEFUL; NUMBER OF INSERTIONS IN THE MIDDLE 80% OF THE GENE (I.E. INSERTIONS IN THE FIRST AND LAST 10% OF THE LENGTH OF THE GENE ARE NOT CONSIDERED)
    sns.violinplot(x='Essentiality',y='Number_Insertions_Truncated_Gene',data=df, cut=0)

    ax2 = sns.boxplot(x='Essentiality',y='Number_Insertions_Truncated_Gene',data=df) #NORMALIZE DATA BY GENE LENGTH
    ax2.set_ylim(-0.001,0.025)
    ax2.set_xlabel('Annotated essential', fontsize=14)
    ax2.set_ylabel('Number of insertions (middle 80% gene) (normalized to gene length)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)


    #NOT USEFUL (?); LARGEST DISTANCE BETWEEN SUBSEQUENT INSERTIONS FOR EACH GENE. Q: WHAT TO DO WITH CASES WHERE THERE IS ONLY A SINGLE OF NO INSERTIONS? -> IF THOSE SITUATIONS SET TO 0 IT DOES GIVE A CLEAR DISTINCTION BETWEEN ESSENTIALITY, BUT IS THIS FAIR?
#    ax = sns.stripplot(x='Essentiality',y='Max_Insertion_Distance', data=df, alpha=0.23, palette='coolwarm')
#
#    sns.violinplot(x='Essentiality',y='Max_Insertion_Distance', data=df, cut=0, palette=['white'])
#
#    df_select = df[df['Number_Insertions_Full_Gene'] > 1]
#    sns.barplot(x='Essentiality',y='Max_Insertion_Distance', data=df_select)
#    del (df_select, ax)


    #POTENTIALLY USEFUL; 
    df_select = df[df['Number_Reads_Full_Gene'] < 10000]
    ax3 = sns.boxplot(x='Essentiality',y='Number_Reads_Full_Gene', data=df_select)
    ax3.set_ylim(-1,800)
    ax3.set_xlabel('Annotated essential', fontsize=14)
    ax3.set_ylabel('Number of reads (normalized to gene length)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    ax4 = sns.boxplot(x='Essentiality',y='Number_Reads_Truncated_Gene', data=df_select)
    ax4.set_ylim(-1,800)
    ax4.set_xlabel('Annotated essential', fontsize=14)
    ax4.set_ylabel('Number of reads (middle 80% gene) (normalized to gene length)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

#    print('Number of outliers for essential genes is %i' % len(boxplot_stats(df.Number_Reads_Truncated_Gene).pop(0)['fliers']))


    ax5 = sns.boxplot(x='Essentiality', y='Reads_per_Transposon_Full_Gene', data=df)
    ax5.set_ylim(-1,100)
    ax5.set_xlabel('Annotated essential', fontsize=14)
    ax5.set_ylabel('Number of reads per transposon (per gene)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)


    ax6 = sns.boxplot(x='Essentiality', y='Reads_per_Transposon_Truncated_Gene', data=df)
    ax6.set_ylim(-1,100)
    ax6.set_xlabel('Annotated essential', fontsize=14)
    ax6.set_ylabel('Number of reads per transposon (per gene, middle 80%)', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    del (df_select, ax1, ax2, ax3, ax4, ax5, ax6)

#%% CREATE SIMILAR GRAPH AS MCIHEL ET.AL. 217 FIG 1E.
#Insertions_NotNorm_Full_Gene
#Reads_NotNorm_Full_Gene
#Number_Insertions_Full_Gene
#Number_Reads_Full_Gene
df_select = df#df[df['Number_Reads_Full_Gene'] < 100]

xmax = 5#100
Nbins = 2000
textsize = 14
line_width = 2
color_nonessential = "#F4A896"
color_essential = "#358597"

fig = plt.figure(figsize=(19,9))
grid = plt.GridSpec(10, 1, wspace=0.0, hspace=0.0)


ax3 = plt.subplot(grid[0,0])
line_props3 = dict(color=color_nonessential, alpha=1.0, linewidth=line_width)
bbox_props3 = dict(color=color_nonessential, alpha=1.0, linewidth=line_width)
median_props3 = dict(color=color_nonessential, linewidth=line_width)
flier_props3 = dict(markeredgecolor=color_nonessential, marker="o", markersize=8)
cap_props3 = dict(color=color_nonessential, linewidth=line_width)
ax3.boxplot(df_select[df_select.Essentiality==False].Number_Insertions_Full_Gene.values, vert=False, whiskerprops=line_props3, boxprops=bbox_props3, flierprops=flier_props3, medianprops=median_props3, capprops=cap_props3, widths=0.9)
ax3.set_xlim(0,xmax)
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.grid(True, linestyle="--", alpha=0.5)


ax4 = plt.subplot(grid[1,0])
line_props4 = dict(color=color_essential, alpha=1.0, linewidth=line_width)
bbox_props4 = dict(color=color_essential, alpha=1.0, linewidth=line_width)
median_props4 = dict(color=color_essential, linewidth=line_width)
flier_props4 = dict(markeredgecolor=color_essential, marker="o", markersize=8)
cap_props4 = dict(color=color_essential, linewidth=line_width)
ax4.boxplot(df_select[df_select.Essentiality==True].Number_Insertions_Full_Gene.values, vert=False, whiskerprops=line_props4, boxprops=bbox_props4, flierprops=flier_props4, medianprops=median_props4, capprops=cap_props4, widths=0.9)
ax4.set_xlim(0,xmax)
ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.grid(True, linestyle="--", alpha=0.5)


ax1 = plt.subplot(grid[2:6,0])
h1, binsize, _ = ax1.hist(df_select[df_select.Essentiality==False].Number_Insertions_Full_Gene.values, bins=Nbins, color=color_nonessential, label="Not annotated essential")
ax1.set_xlim(0,xmax)
ax1.tick_params(labelsize=textsize)
ymax = ax1.get_ylim()
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.set_xticklabels([])


ax2 = plt.subplot(grid[6:11,0])
h2 = ax2.hist(df_select[df_select.Essentiality==True].Number_Insertions_Full_Gene.values, bins=binsize, color=color_essential, label="Annotated essential")
ax2.set_xlim(0,xmax)
ax2.tick_params(labelsize=textsize)
ax2.set_ylim(0,ymax[1])
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.invert_yaxis()


fig.legend(loc="lower right", fontsize=textsize)


#FOLLOWING IS FOR SHOWING COMMON AXIS LABELS.
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("Read density per gene", fontsize=textsize)
plt.ylabel("Number of genes", fontsize=textsize)


plt.tight_layout()


del (fig, grid, ax1, ax2, ax3, ax4, xmax, ymax, Nbins, binsize, textsize, color_nonessential,
     color_essential, line_props3, line_props4, bbox_props3, bbox_props4,
     median_props3, median_props4, flier_props3, flier_props4, cap_props3, cap_props4)


#%% CREATE TEXT FILE FOR NUMBER OF INSERTIONS FOR TRUNCATED GENE
#    savename = "ERR1533148_WT1_insertions_truncated.txt"
#    
#    with open(os.path.join(filepath, savename), 'w') as f:
#        f.write('Gene_name\tEssentiality\tInsertions_10%_truncated\n')
#
#        for gene in gene_inserts_trunc_dict:
#            if gene in essential_position_dict:
#               f. write(gene + '\t' + 'True' + '\t' + str(gene_inserts_trunc_dict.get(gene)).strip('[]') + '\n')
#            else:
#                f.write(gene + '\t' + 'False' + '\t' + str(gene_inserts_trunc_dict.get(gene)).strip('[]') + '\n')
#
#
#    del (savename, f, gene)
#

#%% CREATE TEXT FILE FOR NUMBER OF READS FOR TRUNCATED GENE
#    savename = "ERR1533148_WT1_reads_truncated.txt"
#    
#    with open(os.path.join(filepath, savename), 'w') as f:
#        f.write('Gene_name\tEssentiality\tReads_10%_truncated\n')
#        
#        for gene in gene_reads_trunc_dict:
#                    if gene in essential_position_dict:
#                       f. write(gene + '\t' + 'True' + '\t' + str(gene_reads_trunc_dict.get(gene)).strip('[]') + '\n')
#                    else:
#                        f.write(gene + '\t' + 'False' + '\t' + str(gene_reads_trunc_dict.get(gene)).strip('[]') + '\n')
#
#    del (savename, f, gene)


#%% CREATE SCATTERPLOT NUMBER OF READS PER GENE



#%%
if __name__ == '__main__':
    tninserts_analysis()









