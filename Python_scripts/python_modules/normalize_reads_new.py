# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 09:46:10 2020

@author: gregoryvanbeek
"""
#%%
import os
import numpy as np
file_dirname = os.path.dirname(os.path.abspath('__file__'))
from mapped_reads import total_mapped_reads


#%% TEST PARAMETERS
len_chr = 230218
wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"
#dna_df2; load from genomicfeatures_dataframe_with_normalization.py with parameter 'region=1'


#%%
def reads_normalization_dynamic_window(dna_df2, len_chr, wig_file):
    '''
    '''


    nc_df = dna_df2[dna_df2['Feature_name'] == 'noncoding']
    nc_startposition_list = [position[0] for position in nc_df['Position'].tolist()]
    nc_endposition_list = [position[1] for position in nc_df['Position'].tolist()]
    nc_insrt_list = nc_df['Ninsertions'].tolist()


    for dna in dna_df2.itertuples(index=True):
        if dna.Index < 8:#!!!REMOVE. THIS IS TO STOP AT THE NTH ELEMENT IN DNA_DF2 FOR DEBUGGIN PURPOSES
            print(dna.Feature_name)
            if not dna.Feature_name == 'noncoding':
                
                #GET CLOSEST NONCODING REGION OF CURRENT FEATURE
                feature_start_position = dna.Position[0]
                difference_startposition_featureandnc_list = [abs(nc_startposition - feature_start_position) for nc_startposition in nc_startposition_list]
                index_closest_nc_region = difference_startposition_featureandnc_list.index(min(difference_startposition_featureandnc_list))
                
                #DETERMINE THE NONCODING REGIONS CLOSEST TO THE CURRENT FEATURE
                tn_innc_forwardcount = 0
                index_count = index_closest_nc_region
                for nc_endposition in nc_endposition_list[index_closest_nc_region:]:
                    if tn_innc_forwardcount < 10:
                        print('Forward count nc start position: ', nc_endposition_list[index_count])
                        tn_innc_forwardcount += nc_insrt_list[index_count]
                        index_count += 1
                    else:
                        break
                tn_innc_reversecount = 0
                index_count = 0
                for nc_startposition in nc_startposition_list[:index_closest_nc_region+1]:
                    if tn_innc_reversecount < 10:
                        print('Reverse count nc start position: ', nc_startposition_list[index_count])
                        tn_innc_reversecount += nc_insrt_list[index_count]
                        index_count += 1
                    else:
                        break

        else:#!!!REMOVE
            break#!!!REMOVE
