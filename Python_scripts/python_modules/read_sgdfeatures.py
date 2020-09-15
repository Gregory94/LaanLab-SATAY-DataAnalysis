# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:32:52 2020

@author: gregoryvanbeek

This file reads the SGD_features.txt file found at http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/
"""
#%%
def sgd_features(filepath=None):
    '''
    Output dictionaries contain the following information in this order:
        key:
            0. Feature name
        value:
            0. feature type (l[1])
            1. feature qualifier (Verified or Dubious) (l[2])
            2. Standard name (l[4])
            3. Aliases (separated by '|') (l[5])
            4. Parent feature name (typically 'chromosome ...') (l[6])
            5. Chromosome (l[8])
            6. start coordinate (starting at 0 for each chromosome) (l[9])
            7. end coordinate (starting at 0 for each chromosome) (l[10])
    '''


    if filepath == None:
        filepath = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\SGD_features.tab"


    arabic_to_roman_dict = {1:'I', 2:'II', 3:'III', 4:'IV', 5:'V', 6:'VI',
                            7:'VII', 8:'VIII', 9:'IX', 10:'X', 11:'XI',
                            12:'XII', 13:'XIII', 14:'XIV', 15:'XV', 16:'XVI'}


    with open(filepath) as f:
        lines = f.readlines()


    feature_list = []
    feature_orf_dict = {}
    feature_ars_dict = {}
    feature_telomere_dict = {}
    feature_centromere_dict = {}
    feature_Xelement_dict = {}
    feature_ncrna_dict = {}
    feature_ets_dict = {}
    feature_its_dict = {}
    for line in lines:
        l = line.strip('\n').split('\t')
        if not l[1] in feature_list:
            feature_list.append(l[1])

        if not l[8].endswith('micron') and not l[8] == '':
            chromosome = arabic_to_roman_dict.get(int(l[8]))
            if l[1] == 'ORF':
                feature_orf_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'ARS':
                feature_ars_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'telomere':
                feature_telomere_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'centromere':
                feature_centromere_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'X_element':
                feature_Xelement_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'ncRNA_gene':
                feature_ncrna_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'external_transcribed_spacer_region':
                feature_ets_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]
            elif l[1] == 'internal_transcribed_spacer_region':
                feature_its_dict[l[3]] = [l[1], l[2], l[4], l[5], l[6], chromosome, l[9],l[10]]

    del (lines, f, arabic_to_roman_dict, line, l, chromosome)

    return(feature_orf_dict, feature_ars_dict, feature_telomere_dict,
           feature_centromere_dict, feature_Xelement_dict, feature_ncrna_dict,
           feature_ets_dict, feature_its_dict)

#%%
if __name__ == '_main__':
    sgd_features()
