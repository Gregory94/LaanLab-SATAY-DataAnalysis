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
#ALTERNATIVE FOR CURRENT METHOD: LOOP OVER DNA_DF2, IF FEATURE != NONCODING, GET INDEX IN DNA_DF2 AND FROM THERE ON SEARCH FURTHER FOR ALL NONCODING REGIONS AND REMEMBER THE END POSITION OF THE NONCODING REGION WHEN MORE THAN 10 TRANSPOSONS ARE ENCOUNTERED.
    #DO THE SAME IN REVERSE ORDER TO GET THE PRECEEDING NONCODING REGIONS.

    nc_df = dna_df2[dna_df2['Feature_name'] == 'noncoding']
    nc_startposition_list = [position[0] for position in nc_df['Position'].tolist()]
    nc_endposition_list = [position[1] for position in nc_df['Position'].tolist()]
    nc_insrt_list = nc_df['Ninsertions'].tolist()
    del(nc_df)


    window_start_end_dict = {}
    for dna in dna_df2.itertuples(index=True): # dna = list(dna_df2.itertuples(index=True))[1]
        if not dna.Feature_name == 'noncoding':
            
            #GET CLOSEST NONCODING REGION OF CURRENT FEATURE
            feature_start_position = dna.Position[0]
            difference_startposition_featureandnc_list = [(nc_startposition - feature_start_position) for nc_startposition in nc_startposition_list]
            index_closest_nc_region = np.where(np.array(difference_startposition_featureandnc_list) > 0, difference_startposition_featureandnc_list, np.inf).argmin()
#            index_closest_nc_region = difference_startposition_featureandnc_list.index(min(difference_startposition_featureandnc_list))

            #DETERMINE THE NONCODING REGIONS CLOSEST TO THE CURRENT FEATURE
            tn_innc_forwardcount = 0
            index_count = index_closest_nc_region
            for nc_endposition in nc_endposition_list[index_closest_nc_region:]:
                if tn_innc_forwardcount > 10:
                    window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count-1]]
                    break
                elif index_count >= len(nc_startposition_list)-1:
                    window_start_end_dict[dna.Feature_name] = [nc_endposition_list[index_count]]
                    break
                elif tn_innc_forwardcount < 10:
                    tn_innc_forwardcount += nc_insrt_list[index_count]
                    index_count += 1


#            nc_startposition_list.reverse()
#            tn_innc_reversecount = 0
#            index_count = 0
#            for nc_startposition in nc_startposition_list[index_closest_nc_region-1:]:
#                if tn_innc_reversecount < 10:
#                    print('Reverse count nc start position: ', nc_startposition_list[index_count])
#                    tn_innc_reversecount += nc_insrt_list[index_count]
#                    index_count += 1
#                else:
#                    window_start_end_dict[dna.Feature_name].append(nc_startposition_list[index_count])
#                    break
#                nc_startposition_list.reverse()




