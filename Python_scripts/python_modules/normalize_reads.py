# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:46:10 2020

@author: gregoryvanbeek
"""
#%%
import os
import numpy as np
file_dirname = os.path.dirname(os.path.abspath('__file__'))
#sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from mapped_reads import total_mapped_reads


#%% TEST PARAMETERS
len_chr = 813184#230218
normalization_window_size= 20000
wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
#dna_df2; load from genomicfeatures_dataframe_with_normalization.py with parameter 'region=1'

#%%
def reads_normalization_fixed_window(dna_df2, len_chr, normalization_window_size, wig_file):
    '''
    This function normalizes the number of reads based on windows of fixed size as determined in the normalization_window_size variable.
    It accepts a dataframe that is created using the function dna_features in genomicfeatures_dataframe_with_normalization.py.
    The normalization procedure is based on the method used in Gallagher et.al. 2011.
    This normalizes the reads for the gene length, differences in the number of mapped reads between experiments and the inhomogeniety in the number of insertions within a
    chromosome.
    
    The output of the function is the same dataframe, dna_df2, to which now three more columns are added considering the normalization.
    The first is the normalization described before, the second column is the same normalization, but only for the central part of the gene (for the other genomic features the
    entire feature is taken in consideration) and the third column contains an extra normaliztion where the average number of reads in the noncoding regions in a window are
    defined as fitness 1.
    Also it returns a list including the basepairs where the windows end.
    '''
    

    N = round(len_chr/normalization_window_size)
    window_edge_list = np.linspace(0, len_chr, N, dtype=int).tolist()#[82500, 243500, len_chr]
    window_length = window_edge_list[1] - window_edge_list[0] + 1

    total_reads_in_genome = total_mapped_reads(wig_file)


#    read_density_chromosome = sum([(dna.Nreads/dna.Nbasepairs) for dna in dna_df2.itertuples() if dna.Feature_type == None])
    read_density_chromosome = sum([dna.Nreads for dna in dna_df2.itertuples() if dna.Feature_type == None]) / sum([dna.Nbasepairs for dna in dna_df2.itertuples() if dna.Feature_type == None])
    read_density_windows = np.ones(len(window_edge_list)-1)
    window_start = 0
    i = 0
    for window_end in window_edge_list[1:]:
        read_density = sum([(dna.Nreads) for dna in dna_df2.itertuples() if window_start <= dna.Position[0] < window_end and dna.Feature_type == None]) / sum([(dna.Nbasepairs) for dna in dna_df2.itertuples() if window_start <= dna.Position[0] < window_end and dna.Feature_type == None])
        if not read_density == 0:
            read_density_windows[i] = read_density
#            read_density_windows.append(sum([nc[14] for nc in dna_df2.itertuples() if window_start <= nc[6][0] < window_end and nc[4] == None]))
        window_start = window_end
        i += 1


    norm_reads_list = []
    norm_reads_truncatedgene_list = []
    window_index_list = [] #LIST CONTAINING THE INDICES FOR EACH FEATURE INDICATING TO WHICH WINDOW IT BELONGS. THIS HAS THE SAME LENGTH AS DNA_DF2 WHERE EACH ROW CORRESPONDS TO THE ROW IN DNA_DF2
    for row in dna_df2.itertuples():
        #normalization equation:
            #normalization for gene = raw read count in middle 80% of feature * (1/gene length) * (10^6/total mapped reads in genome) * (read density in chromosome/read density in window)
            #normalization for other features = raw read count in entire feature * (1/gene length) * (10^6/total mapped reads in genome) * (read density in chromosome/read density in window)
        read_density_windows_index = int(row.Position[0]/window_length) #determine which window the current feature belongs to.
        window_index_list.append(read_density_windows_index)
        if not row.Feature_type == None and row.Feature_type.startswith('Gene'):
            norm_reads_truncatedgene_list.append(row.Nreads_truncatedgene * (1/row.Nbasepairs*0.8) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/read_density_windows[read_density_windows_index]))
        else:
            norm_reads_truncatedgene_list.append(row.Nreads_truncatedgene * (1/row.Nbasepairs*1.0) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/read_density_windows[read_density_windows_index]))
        norm_reads_list.append(row.Nreads * (1/row.Nbasepairs) * ((10**6)/(total_reads_in_genome)*1.0) * (read_density_chromosome/read_density_windows[read_density_windows_index]))
        

    dna_df2['Nreads_normalized'] = norm_reads_list
    dna_df2['Nreads_truncatedgene_normalized'] = norm_reads_truncatedgene_list


    #NORMALIZE EACH WINDOW WITH THE AVERAGE NUMBER OF READS IN THE NONCODING REGIONS.
    readsperbp_noncoding = [0] * (len(window_edge_list) - 1) #CONTAINS THE SUM OF THE NUMBER OF READS IN THE NONCODING REGIONS FOR EACH WINDOW.
    N_reads_noncoding = [0] * (len(window_edge_list) - 1) #CONTAINS THE NUMBER OF NONCODING REGIONS IN EACH WINDOW.
    i = 0
    for row in dna_df2.itertuples():
        if row.Feature_type == None:
            readsperbp_noncoding[window_index_list[i]] += row.Nreads_normalized
            N_reads_noncoding[window_index_list[i]] += 1
        i += 1
    avg_reads_noncoding = [readsperbp_noncoding[i] / N_reads_noncoding[i] for i in range(len(readsperbp_noncoding))]
    
    normbync_reads_list = []
    for row in dna_df2.itertuples():
        normbync_reads_list.append(row.Nreads_normalized / avg_reads_noncoding[int(row.Position[0] / window_length)])
    dna_df2['Nreads_normalized_byNCregions'] = normbync_reads_list



    return(dna_df2, window_edge_list)


#%%
def reads_normalization_dynamic_window(dna_df2, len_chr, normalization_window_size, wig_file):
    '''
    This function normalizes the number of reads based on windows of dynamic size.
    Each window has at most a number of basepairs as determined in normalization_window_size or a maximum of about 10 transposon insertions.
    It accepts a dataframe that is created using the function dna_features in genomicfeatures_dataframe_with_normalization.py.
    The normalization procedure is based on the method used in Gallagher et.al. 2011.
    This normalizes the reads for the gene length, differences in the number of mapped reads between experiments and the inhomogeniety in the number of insertions within a
    chromosome.
    
    The output of the function is the same dataframe, dna_df2, to which now three more columns are added considering the normalization.
    The first is the normalization described before, the second column is the same normalization, but only for the central part of the gene (for the other genomic features the
    entire feature is taken in consideration) and the third column contains an extra normaliztion where the average number of reads in the noncoding regions in a window are
    defined as fitness 1.
    Also it returns a list including the basepairs where the windows end.
    '''

#DETERMINE THE EDGES OF THE WINDOWS
#!!! CREATE NEW CODE THAT TAKES THE GENES, AND NORMALIZE WITH NC REGIONS BEFORE AND AFTER THAT ADDS UP TO N INSERTIONS -> SEE NORMALIZE_READS_NEW.PY
#    bp_count = 0
    tn_count = 0
    window_edge_list = []
    for dna in dna_df2.itertuples():
        if dna.Feature_type == None:
#            bp_count += dna.Nbasepairs
            tn_count += dna.Ninsertions
        if tn_count > 30:
            window_edge_list.append(dna.Position[1])
#            bp_count = 0
            tn_count = 0
    window_edge_list[-1] = dna_df2.iloc[-1].Position[1]


#FOR EACH WINDOW, DETERMINE ITS SIZE
    window_size_list = [window_edge_list[0]]
    for i in range(len(window_edge_list[1:])):
        window_size_list.append(window_edge_list[i+1] - window_edge_list[i])

#GET THE TOTAL NUMBER OF READS IN THE GENOME
    total_reads_in_genome = total_mapped_reads(wig_file)


#
    read_density_chromosome = sum([dna.Nreads for dna in dna_df2.itertuples() if dna.Feature_type == None]) / sum([dna.Nbasepairs for dna in dna_df2.itertuples() if dna.Feature_type == None])
    read_density_windows = np.ones(len(window_edge_list))
    window_start = 0
    i = 0
    for window_end in window_edge_list:
        read_density = sum([(dna.Nreads) for dna in dna_df2.itertuples() if window_start <= dna.Position[0] < window_end and dna.Feature_type == None]) / sum([(dna.Nbasepairs) for dna in dna_df2.itertuples() if window_start <= dna.Position[0] < window_end and dna.Feature_type == None])
        if not read_density == 0:
            read_density_windows[i] = read_density
        window_start = window_end
        i += 1


    norm_reads_list = []
    norm_reads_truncatedgene_list = []
    window_index_list = [] #LIST CONTAINING THE INDICES FOR EACH FEATURE INDICATING TO WHICH WINDOW IT BELONGS. THIS HAS THE SAME LENGTH AS DNA_DF2 WHERE EACH ROW CORRESPONDS TO THE ROW IN DNA_DF2
    for row in dna_df2.itertuples():
        for ind, w in reversed(list(enumerate(window_edge_list))):
            if row.Position[0] < window_edge_list[0]:
                window_length = window_size_list[0]
                ind = -1
                break
            elif row.Position[0] > w:
                window_length = window_size_list[ind+1]
                break
        
        read_density_windows_index = ind + 1
        window_index_list.append(read_density_windows_index)
        if not row.Feature_type == None and row.Feature_type.startswith('Gene'):
            norm_reads_truncatedgene_list.append(row.Nreads_truncatedgene * (1/row.Nbasepairs*0.8) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/read_density_windows[read_density_windows_index]))
        else:
            norm_reads_truncatedgene_list.append(row.Nreads_truncatedgene * (1/row.Nbasepairs*1.0) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/read_density_windows[read_density_windows_index]))
        norm_reads_list.append(row.Nreads * (1/row.Nbasepairs) * ((10**6)/(total_reads_in_genome)*1.0) * (read_density_chromosome/read_density_windows[read_density_windows_index]))


    dna_df2['Nreads_normalized'] = norm_reads_list
    dna_df2['Nreads_truncatedgene_normalized'] = norm_reads_truncatedgene_list


    #NORMALIZE EACH WINDOW WITH THE AVERAGE NUMBER OF READS IN THE NONCODING REGIONS.
    readsperbp_noncoding = [0] * (len(window_edge_list)) #CONTAINS THE SUM OF THE NUMBER OF READS IN THE NONCODING REGIONS FOR EACH WINDOW.
    N_reads_noncoding = [0] * (len(window_edge_list)) #CONTAINS THE NUMBER OF NONCODING REGIONS IN EACH WINDOW.
    i = 0
    for row in dna_df2.itertuples():
        
        if row.Feature_type == None:
            readsperbp_noncoding[window_index_list[i]] += row.Nreads_normalized
            N_reads_noncoding[window_index_list[i]] += 1
        i += 1
    avg_reads_noncoding = [readsperbp_noncoding[i] / N_reads_noncoding[i] for i in range(len(readsperbp_noncoding))]


    normbync_reads_list = []
    for row in dna_df2.itertuples():
        for ind, w in reversed(list(enumerate(window_edge_list))):
            if row.Position[0] < window_edge_list[0]:
                window_length = window_size_list[0]
                ind = -1
                break
            elif row.Position[0] > w:
                window_length = window_size_list[ind+1]
                break
        normbync_reads_list.append(row.Nreads_normalized / avg_reads_noncoding[ind+1])#int(row.Position[0] / window_length)])
    dna_df2['Nreads_normalized_byNCregions'] = normbync_reads_list



    return(dna_df2, window_edge_list)