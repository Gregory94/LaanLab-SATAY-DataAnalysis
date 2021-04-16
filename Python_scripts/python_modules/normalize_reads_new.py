# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:46:10 2020

@author: gregoryvanbeek
"""
#%%
#import os
import numpy as np
#file_dirname = os.path.dirname(os.path.abspath('__file__'))
from mapped_reads import total_mapped_reads


#%% TEST PARAMETERS
len_chr = 316620#230218
wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
Ninsrt_threshold = 20 #NUMBER OF INSERTIONS BEFORE AND AFTER A FEATURE. THIS IS NOT THE NUMBER OF READS.
#dna_df2; load from genomicfeatures_dataframe_with_normalization.py with parameter 'region=1'


#%%
def reads_normalization_dynamic_window(dna_df2, len_chr, wig_file, Ninsrt_threshold):
    '''
    This function normalizes the read count in all features in a chromosome.
    It loops over all features in dna_df2 (created in genomicfeatures_dataframe.py).
    For each feature it looks for all flanking regions upstream and downstream until the number of insertions exceed Ninsrt_threshold.
    All flanking regions define the window of a specific feature and all the noncoding regions in the window are used to determine the read density.
    This read density is used for normalization in the following equation:
        norm_read = sum all read in feature * (1/gene length) * (1/total mapped reads in genome) * (read density chromosome / read density window)
    The output is added as a column in dna_df2 with the header 'Nreads_truncatedgene_normalized'.
    '''
    nc_df = dna_df2[dna_df2['Feature_name'] == 'noncoding']
    nc_startposition_list = [position[0] for position in nc_df['Position'].tolist()]
    nc_endposition_list = [position[1] for position in nc_df['Position'].tolist()]
    nc_insrt_list = nc_df['Ninsertions'].tolist()
    nc_read_list = nc_df['Nreads'].tolist()
    del(nc_df)


    total_reads_in_genome = total_mapped_reads(wig_file)
    read_density_chromosome = sum([dna.Nreads for dna in dna_df2.itertuples() if dna.Feature_type == None]) / sum([dna.Nbasepairs for dna in dna_df2.itertuples() if dna.Feature_type == None])


    #IF LENGTH OF DICTIONARIES IS DIFFERENT THAN LENGTH DNA_DF2, IT MEANS THAT SOME FEATURE APPEAR MULTIPLE TIMES IN DNA_DF2. THE KEYS IN THE DICTIONARIES ARE THEN OVERWRITTEN.
    window_start_end_dict = {} #-> EACH VALUE CONTAINS TWO INTEGERS IN A LIST: 1;START_POSITION NONCODING REGION BEFORE CURRENT GENE THAT HAS MORE THAN N TRANSPOSONS 2;END_POSITION NONCODING REGION AFTER CURRENT GENE THAT HAS MORE THAN N TRANPOSONS
    window_insrt_dict = {} #-> KEYS ARE THE SAME AS WINDOW_START_END_DICT AND VALUES CONTAIN TWO INTEGERS: 1;NUMBER OF INSERTIONS DOWNSTREAM OF GENE 2;NUMBER OF INSERTIONS UPSTREAM OF GENE.
    window_read_dict = {} #-> KEYS ARE THE SAME AS WINDOW_START_END_DICT AND VALUES CONTAIN TOTAL NUMBER OF READS DOWNSTREAM AND UPSTREAM OF GENE.
    window_nclength_dict = {} #-> KEYS ARE THE SAME AS WINDOW_START_END_DICT AND VALUE CONTAIN THE TOTAL LENGTH OF ALL NONCODING REGIONS IN THE WINDOW
    for dna in dna_df2.itertuples(index=True): # dna = list(dna_df2.itertuples(index=True))[1]

        #GET CLOSEST NONCODING REGION OF CURRENT FEATURE
        feature_start_position = dna.Position[0]
        difference_startposition_featureandnc_list = [(nc_startposition - feature_start_position) for nc_startposition in nc_startposition_list]
        index_closest_nc_region = np.where(np.array(difference_startposition_featureandnc_list) > 0, difference_startposition_featureandnc_list, np.inf).argmin()

        #DETERMINE THE NONCODING REGIONS CLOSEST TO THE CURRENT FEATURE (reversecount <- gene -> forwardcount)
        tn_innc_bpcount = 0
        read_innc_count = 0
        if index_closest_nc_region == 0 and not dna.Feature_name == 'noncoding':
            window_start_end_dict[dna.Feature_name] = [len_chr]
            window_insrt_dict[dna.Feature_name] = [0]
        elif index_closest_nc_region == 0 and dna.Feature_name == 'noncoding':
            window_start_end_dict[dna.Feature_name + str(dna.Position[0])] = [len_chr]
            window_insrt_dict[dna.Feature_name + str(dna.Position[0])] = [0]
        elif not index_closest_nc_region == 0:
            tn_innc_forwardcount = 0
            index_count = index_closest_nc_region
            for nc_endposition in nc_endposition_list[index_closest_nc_region:]:
                if tn_innc_forwardcount > Ninsrt_threshold:
                    if not dna.Feature_name == 'noncoding':
                        window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count-1]]
                        window_insrt_dict[dna.Feature_name] = [tn_innc_forwardcount]
                    elif dna.Feature_name == 'noncoding':
                        window_start_end_dict[dna.Feature_name + str(dna.Position[0])] = [nc_endposition_list[index_count-1]]
                        window_insrt_dict[dna.Feature_name + str(dna.Position[0])] = [tn_innc_forwardcount]
                    break
                elif index_count >= len(nc_endposition_list)-1:
                    if not dna.Feature_name == 'noncoding':
                        window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count]]
                        window_insrt_dict[dna.Feature_name] = [tn_innc_forwardcount]
                    elif dna.Feature_name == 'noncoding':
                        window_start_end_dict[dna.Feature_name + str(dna.Position[0])] = [nc_endposition_list[index_count]]
                        window_insrt_dict[dna.Feature_name + str(dna.Position[0])] = [tn_innc_forwardcount]
                    break
                elif tn_innc_forwardcount <= Ninsrt_threshold:
                    tn_innc_forwardcount += nc_insrt_list[index_count]
                    tn_innc_bpcount += nc_endposition_list[index_count] - nc_startposition_list[index_count] + 1
                    read_innc_count += nc_read_list[index_count]
                    index_count += 1

        tn_innc_reversecount = 0
        if index_closest_nc_region == 0:
            index_closest_nc_region = len(nc_read_list)-1
            index_count = index_closest_nc_region
        elif not index_closest_nc_region == 0:
            index_count = index_closest_nc_region-1
        for nc_startposition in nc_startposition_list[:index_closest_nc_region]:
            if tn_innc_reversecount > Ninsrt_threshold:
                if not dna.Feature_name == 'noncoding':
                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count+1]] + window_start_end_dict.get(dna.Feature_name)
                    window_insrt_dict[dna.Feature_name] = [tn_innc_reversecount] + window_insrt_dict.get(dna.Feature_name)
                    window_nclength_dict[dna.Feature_name] = tn_innc_bpcount
                    window_read_dict[dna.Feature_name] = read_innc_count
                elif dna.Feature_name == 'noncoding':
                    window_start_end_dict[dna.Feature_name + str(dna.Position[0])] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name + str(dna.Position[0]))
                    window_insrt_dict[dna.Feature_name + str(dna.Position[0])] = [tn_innc_reversecount] + window_insrt_dict.get(dna.Feature_name + str(dna.Position[0]))
                    window_nclength_dict[dna.Feature_name + str(dna.Position[0])] = tn_innc_bpcount
                    window_read_dict[dna.Feature_name + str(dna.Position[0])] = read_innc_count
                break
            elif index_count <= 1:
                if not dna.Feature_name == 'noncoding':
                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name)
                    window_insrt_dict[dna.Feature_name] = [tn_innc_reversecount] + window_insrt_dict.get(dna.Feature_name)
                    window_nclength_dict[dna.Feature_name] = tn_innc_bpcount
                    window_read_dict[dna.Feature_name] = read_innc_count
                elif dna.Feature_name == 'noncoding':
                    window_start_end_dict[dna.Feature_name + str(dna.Position[0])] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name + str(dna.Position[0]))
                    window_insrt_dict[dna.Feature_name + str(dna.Position[0])] = [tn_innc_reversecount] + window_insrt_dict.get(dna.Feature_name + str(dna.Position[0]))
                    window_nclength_dict[dna.Feature_name + str(dna.Position[0])] = tn_innc_bpcount
                    window_read_dict[dna.Feature_name + str(dna.Position[0])] = read_innc_count
                break
            elif tn_innc_reversecount <= Ninsrt_threshold:
                if not dna.Feature_name == 'noncoding':
                    tn_innc_reversecount += nc_insrt_list[index_count]
                    tn_innc_bpcount += nc_endposition_list[index_count] - nc_startposition_list[index_count] + 1
                    read_innc_count += nc_read_list[index_count]
                elif dna.Feature_name == 'noncoding':
                    tn_innc_reversecount += nc_insrt_list[index_count-1]
                    tn_innc_bpcount += nc_endposition_list[index_count-1] - nc_startposition_list[index_count-1] + 1
                    read_innc_count += nc_read_list[index_count-1]
                index_count -= 1


    del (difference_startposition_featureandnc_list, dna, feature_start_position, index_closest_nc_region, index_count, nc_endposition, nc_startposition, tn_innc_forwardcount, tn_innc_reversecount, nc_insrt_list, tn_innc_bpcount)



    #!!!HOW TO DEAL WITH THE DIFFERENCE IN LENGTH OF DICTIONARY AND DNA_DF2? -> CHANGE KEYS DICTIONARIES FROM GENE NAMES TO DNA_DF2 INDEX?
    norm_reads_list = []
    for dna in dna_df2.itertuples():
        if not dna.Feature_name == 'noncoding' and dna.Feature_type.startswith('Gene'):
            norm_reads_list.append(dna.Nreads_truncatedgene * (1/dna.Nbasepairs*0.8) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/(window_read_dict.get(dna.Feature_name) / window_nclength_dict.get(dna.Feature_name))))
        elif not dna.Feature_name == 'noncoding':
            norm_reads_list.append(dna.Nreads_truncatedgene * (1/dna.Nbasepairs*1.0) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/(window_read_dict.get(dna.Feature_name) / window_nclength_dict.get(dna.Feature_name))))
        elif dna.Feature_name == 'noncoding':
            norm_reads_list.append(dna.Nreads_truncatedgene * (1/dna.Nbasepairs*1.0) * ((10**6)/(total_reads_in_genome)) * (read_density_chromosome/(window_read_dict.get(dna.Feature_name + str(dna.Position[0])) / window_nclength_dict.get(dna.Feature_name + str(dna.Position[0])))))

    dna_df2['Nreads_truncatedgene_normalized'] = norm_reads_list
    
    return(dna_df2, nc_startposition_list, nc_endposition_list)














