# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:46:10 2020

@author: gregoryvanbeek
"""
#%%
import os
import numpy as np
#file_dirname = os.path.dirname(os.path.abspath('__file__'))
#from mapped_reads import total_mapped_reads


#%% TEST PARAMETERS
len_chr = 230218
wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
#dna_df2; load from genomicfeatures_dataframe_with_normalization.py with parameter 'region=1'


#%%
def reads_normalization_dynamic_window(dna_df2, len_chr, wig_file):
    '''
    '''
    Ninsrt_threshold = 10

    nc_df = dna_df2[dna_df2['Feature_name'] == 'noncoding']
    nc_startposition_list = [position[0] for position in nc_df['Position'].tolist()]
    nc_endposition_list = [position[1] for position in nc_df['Position'].tolist()]
    nc_insrt_list = nc_df['Ninsertions'].tolist()
    del(nc_df)


    window_start_end_dict = {} #-> CONTAINS TWO INTEGERS: 1;START_POSITION NONCODING REGION BEFORE CURRENT GENE THAT HAS MORE THAN N TRANSPOSONS 2;END_POSITION NONCODING REGION AFTER CURRENT GENE THAT HAS MORE THAN N TRANPOSONS
    for dna in dna_df2.itertuples(index=True): # dna = list(dna_df2.itertuples(index=True))[1]
        if not dna.Feature_name == 'noncoding':

            #GET CLOSEST NONCODING REGION OF CURRENT FEATURE
            feature_start_position = dna.Position[0]
            difference_startposition_featureandnc_list = [(nc_startposition - feature_start_position) for nc_startposition in nc_startposition_list]
            index_closest_nc_region = np.where(np.array(difference_startposition_featureandnc_list) > 0, difference_startposition_featureandnc_list, np.inf).argmin()

            #DETERMINE THE NONCODING REGIONS CLOSEST TO THE CURRENT FEATURE (reversecount <- gene -> forwardcount)
            tn_innc_forwardcount = 0
            index_count = index_closest_nc_region
            for nc_endposition in nc_endposition_list[index_closest_nc_region:]:
                if tn_innc_forwardcount > Ninsrt_threshold:
                    window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count-1]]
                    break
                elif index_count >= len(nc_endposition_list)-1:
                    window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count]]
                    break
                elif tn_innc_forwardcount <= Ninsrt_threshold:
                    tn_innc_forwardcount += nc_insrt_list[index_count]
                    index_count += 1

            tn_innc_reversecount = 0
            index_count = index_closest_nc_region
            for nc_startposition in nc_startposition_list[:index_closest_nc_region]:
                if tn_innc_reversecount > Ninsrt_threshold:
                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name)
                    break
                elif index_count <= 1:
                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count-1]] + window_start_end_dict.get(dna.Feature_name)
                    break
                elif tn_innc_reversecount <= Ninsrt_threshold:
                    tn_innc_reversecount += nc_insrt_list[index_count]
                    index_count -= 1


##            if index_closest_nc_region == 1:
##                window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_closest_nc_region-1]] + window_start_end_dict.get(dna.Feature_name)
##            else:
#            tn_innc_reversecount = 0
#            index_count = index_closest_nc_region
#            for inx in range(1, len(nc_startposition_list[:index_closest_nc_region])):
#                if index_count <= 1:
#                    print(dna.Feature_name)
#                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name)
#                    break
#                elif tn_innc_reversecount > Ninsrt_threshold:
#                    window_start_end_dict[dna.Feature_name] = [nc_startposition_list[index_count]] + window_start_end_dict.get(dna.Feature_name)
#                    break
#                elif tn_innc_reversecount <= Ninsrt_threshold:
#                    tn_innc_reversecount += nc_insrt_list[index_count]
#                    index_count -= 1
##                    inx -= 1





