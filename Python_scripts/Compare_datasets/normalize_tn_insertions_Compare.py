# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:02:44 2020

@author: gregoryvanbeek

This file can be found at https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/dev_Gregory/Python_scripts/normalize_tn_insertions_Compare.py

This script is based on 'normalize_tn_insertions.py which can be found at:
    https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/blob/dev_Gregory/Python_scripts/normalize_tn_insertions.py
"""

#%%
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'..','python_modules'))
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic
from chromosome_names_in_files import chromosome_name_wigfile
from gene_names import list_gene_names, gene_aliases
from read_sgdfeatures import sgd_features


#%%
def dna_features_compare(region, wig_file_list, pergene_insertions_file_list):
    '''This function inputs two wig files and two pergene_insertions files created using transposonmapping_satay.py.
    The The two files are compared with each other in terms of number of insertions per basepair per genomic region.
    Output is a barplot indicating the relative difference in number of transposon insertions per basepair per genomic region.
    This is determined as the the number of insertions per basepair in the first file - the number of insertions per basepair in the second file divided by the maximum number of insertions per basepair in the first file.
    A genomic region is here defined as a gene (separated as annotated essential and not essential), telomere, centromere, ars etc.
    The width of the bars correspond to the width of the genomic region.
    
    
    TO USE THIS SCRIPT:
        Change the input file paths in the last lines of this code.
        Change the input files in the next section (to get the files, check https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/dev_Gregory/Python_scripts/Data_Files)
        Make sure that the imported modules (previous section) are all present in a folder called 'python_modules' that is at the same location as where this script is saved (see https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/dev_Gregory/Python_scripts/python_modules).
        
        All information that is used to create the plots in this script are saved in the dataframe variable 'dna_df2'.
    '''
#%% INPUT FILES
    essentials_file = r"C:\Users\linigodelacruz\Documents\PhD_2018\Documentation\SATAY\src(source-code)\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_AllEssentialGenes_List.txt"

    gene_information_file = os.path.join(file_dirname,'..','Data_Files','Yeast_Protein_Names.txt')

    gff_file = os.path.join(file_dirname,'..','Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')

    sgd_features_file = os.path.join(file_dirname,'..','Data_Files','SGD_features.tab')

#%% DETERMINE INPUTTED REGION

    if type(region) == str:
        if region in chromosomename_roman_to_arabic()[1]:
            chrom = region
            region_start = None
            region_end = None
        elif region in list_gene_names(gene_information_file):
            #LOOK UP GENE INFORMATION ...
            pass
        else:
            print("WARNING: Specified chromosome not found. Enter chromosome as a roman numeral between I and XVI")

    elif type(region) == list: #ADD CODE THAT ONLY TAKES GENES WITHIN SPECIFIED REGION
        chrom = region[0]
        region_start = region[1]
        region_end = region[2]


    del (region)

#%% READ WIG FILE FOR GETTING LOCATIONS OF ALL TN INSERTIONS
    for file_counter in range(len(wig_file_list)):
        print("Reading file: %s" % (wig_file_list[file_counter]))
        with open(wig_file_list[file_counter], 'r') as f:
            lines = f.readlines()
    
        chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(lines)[1:]
    
        insrt_in_chrom_list = []
        reads_in_chrom_list = []
        for l in lines[chrom_start_line_dict.get(chrom):chrom_end_line_dict.get(chrom)]:
            insrt_in_chrom_list.append(int(l.strip('\n').split(' ')[0]))
            reads_in_chrom_list.append(int(l.strip('\n').split(' ')[1]))
    
    
        del (lines, l, f, chrom_start_line_dict, chrom_end_line_dict)
    
    #%% READ PERGENE_INSERTIONS FILE FOR LOCATION OF ALL INSERTIONS PER EACH GENE.
    
        with open(pergene_insertions_file_list[file_counter]) as f:
            lines = f.readlines()
    
    
        gene_position_dict = {}
        gene_inserts_dict = {}
        gene_reads_dict = {}
        for line in lines[1:]:
            line_split = line.strip('\n').split('\t')
    
    
            if line_split[1] == chrom:
                genename = line_split[0]
                gene_chrom = line_split[1]
                gene_start = int(line_split[2])
                gene_end = int(line_split[3])
    
                gene_position_dict[genename] = [gene_chrom, gene_start, gene_end] #DICT CONTAINING ALL GENES WITHIN THE DEFINED CHROMOSOME INCLUDING ITS START AND END POSITION
    
    
                geneinserts_str = line_split[4].strip('[]')
                if not geneinserts_str == '':
                    geneinserts_list = [int(ins) for ins in geneinserts_str.split(',')]
                else:
                    geneinserts_list = []
    
                gene_inserts_dict[genename] = geneinserts_list
    
    
                genereads_str = line_split[5].strip('[]')
                if not genereads_str == '':
                    genereads_list = [int(read) for read in genereads_str.split(',')]
                else:
                    genereads_list = []
    
                gene_reads_dict[genename] = genereads_list
    
    
                if len(geneinserts_list) != len(genereads_list):
                    print('WARNING: %s has different number of reads compared with the number of inserts' % genename )
    
    
        del (f, lines, line, genename, gene_chrom, gene_start, gene_end, geneinserts_list, geneinserts_str, genereads_str, genereads_list)
    
    #%% DETERMINE THE LOCATION GENOMIC FEATURES IN THE CURRENT CHROMOSOME AND STORE THIS IN A DICTIONARY
        len_chr = chromosome_position(gff_file)[0].get(chrom)
        start_chr = chromosome_position(gff_file)[1].get(chrom)
        end_chr = chromosome_position(gff_file)[2].get(chrom)
    
        dna_dict = {} #for each bp in chromosome, determine whether it belongs to a noncoding or coding region
        for bp in range(start_chr, end_chr + 1): #initialize dna_dict
            dna_dict[bp] = ['noncoding', None] #form is: ['element_name', 'type']
    
    
    
        feature_orf_dict = sgd_features(sgd_features_file)[0]
        gene_alias_dict = gene_aliases(gene_information_file)[0]
    
        for gene in gene_position_dict:
            if gene in feature_orf_dict:
                for bp in range(gene_position_dict.get(gene)[1]+start_chr, gene_position_dict.get(gene)[2]+start_chr+1):
                    dna_dict[bp] = [gene, "Gene; "+feature_orf_dict.get(gene)[1]]
            else:
                gene_alias = [key for key, val in gene_alias_dict.items() if gene in val][0]
                for bp in range(gene_position_dict.get(gene)[1]+start_chr, gene_position_dict.get(gene)[2]+start_chr+1):
                    dna_dict[bp] = [gene_alias, "Gene; "+feature_orf_dict.get(gene_alias)[1]]
    
        del (gene, bp)
    
    #%% GET FEATURES FROM INTERGENIC REGIONS (-> SEE SGD_features.tab IN DATA_FILES IN GITHUB FOLDER) -> MAKE FUNCTION FOR THIS
    
        dna_dict = feature_position(sgd_features(sgd_features_file)[1], chrom, start_chr, dna_dict, "ARS") #ARS
        dna_dict = feature_position(sgd_features(sgd_features_file)[2], chrom, start_chr, dna_dict, "Telomere") #Telomeres
        dna_dict = feature_position(sgd_features(sgd_features_file)[3], chrom, start_chr, dna_dict, "Centromere") #Centromeres
        dna_dict = feature_position(sgd_features(sgd_features_file)[4], chrom, start_chr, dna_dict, "X_element") #X_element
        dna_dict = feature_position(sgd_features(sgd_features_file)[5], chrom, start_chr, dna_dict, "ncRNA") #ncRNA
        dna_dict = feature_position(sgd_features(sgd_features_file)[6], chrom, start_chr, dna_dict, "External_transcribed_Spacer_Region") #External transcribed spacer region
        dna_dict = feature_position(sgd_features(sgd_features_file)[7], chrom, start_chr, dna_dict, "Internal_transcribed_Spacer_Region") #Internal transcribed spacer region
    
    
        ### TEST IF ELEMENTS IN FEATURE_ORF_DICT FOR SELECTED CHROMOSOME ARE THE SAME AS THE GENES IN GENE_POSITION_DICT BY CREATING THE DICTIONARY FEATURE_POSITION_DICT CONTAINING ALL THE GENES IN FEATURE_ORF_DICT WITH THEIR CORRESPONDING POSITION IN THE CHROMOSOME
        gene_alias_dict = gene_aliases(gene_information_file)[0]
        orf_position_dict = {}
        for feature in feature_orf_dict:
            if feature_orf_dict.get(feature)[5] == chrom:
                if feature in gene_position_dict:
                    orf_position_dict[feature] = [feature_orf_dict.get(feature)[6], feature_orf_dict.get(feature)[7]]
                else:
                    for feature_alias in gene_alias_dict.get(feature):
                        if feature_alias in gene_position_dict:
                            orf_position_dict[feature_alias] = [feature_orf_dict.get(feature)[6], feature_orf_dict.get(feature)[7]]
    
    
    
        if sorted(orf_position_dict) == sorted(gene_position_dict):
            print('Everything alright, just ignore me!')
        else:
            print('WARNING: Genes in feature_list are not the same as the genes in the gene_position_dict. Please check!')
    
    
        del (feature_orf_dict, orf_position_dict, feature, feature_alias)
    
    #%% DETERMINE THE NUMBER OF TRANSPOSONS PER BP FOR EACH FEATURE AND STORE RESULTS IN DATAFRAME
    
    
        dna_df = pd.DataFrame.from_dict(dna_dict, orient="index")
    
    
        reads_loc_list = [0] * len(dna_dict) # CONTAINS ALL READS JUST LIKE READS_IN_CHROM_LIST, BUT THIS LIST HAS THE SAME LENGTH AS THE NUMBER OF BP IN THE CHROMOSOME WHERE THE LOCATIONS WITH NO READS ARE FILLED WITH ZEROS
        i = 0
        for ins in insrt_in_chrom_list:
            reads_loc_list[ins] = reads_in_chrom_list[i]
            i += 1
        dna_df["Reads"] = reads_loc_list
    
    
        del (i, ins, dna_df)
    
    #%% CREATE DATAFRAME FOR EACH FEATURE (E.G. NONCODING DNA, GENE, ETC.) IN THE CHROMOSOME AND DETERMINE THE NUMBER OF INSERTIONS AND READS PER FEATURE.
        feature_NameAndType_list = []
        f_previous = dna_dict.get(start_chr)[0]
        f_type = dna_dict.get(start_chr)[1]
        N_reads = 0
        N_reads_list = []
        N_insrt = 0
        N_insrt_list = []
        N_bp = 1
        N_bp_list = []
        i = 0
        for bp in dna_dict:
            f_current = dna_dict.get(bp)[0]
            if f_current == f_previous:
                f_type = dna_dict.get(bp)[1]
                N_bp += 1
                N_reads += reads_loc_list[i]
                if not reads_loc_list[i] == 0:
                    N_insrt += 1
            elif (f_current != f_previous or (i+start_chr) == end_chr):# and not f_current.endswith('-A'):
                feature_NameAndType_list.append([f_previous, f_type])
                N_reads_list.append(N_reads)
                N_insrt_list.append(N_insrt)
                N_bp_list.append(N_bp)
                N_reads = 0
                N_insrt = 0
                N_bp = 1
                f_previous = f_current
            i += 1
    
        N_reads_per_bp_list = []
        N_insrt_per_bp_list = []
        for i in range(len(N_reads_list)):
            N_reads_per_bp_list.append(N_reads_list[i]/N_bp_list[i])
            N_insrt_per_bp_list.append(N_insrt_list[i]/N_bp_list[i])
    
    
        #############get all essential genes together with their aliases##############
        with open(essentials_file, 'r') as f:
            essentials_temp_list = f.readlines()[1:]
        essentials_list = [essential.strip('\n') for essential in essentials_temp_list]
        del essentials_temp_list
    
        gene_alias_dict = gene_aliases(gene_information_file)[0]
        for key, val in gene_alias_dict.items():
            if key in essentials_list:
                for alias in val:
                    essentials_list.append(alias)
    
        #ADD
        essentiality_list = []
        for feature in feature_NameAndType_list:
            if not feature[0] == "noncoding":
                if feature[1] == "ARS" or feature[1] == "Telomere" or feature[1] == "Centromere" or feature[1] == "X_element" or feature[1] == "ncRNA" or feature[1] == "External_transcribed_Spacer_Region" or feature[1] == "Internal_transcribed_Spacer_Region":
                    essentiality_list.append(None)
                elif feature[0] in essentials_list:
                    essentiality_list.append(True)
                else:
                    essentiality_list.append(False)
            else:
                essentiality_list.append(None)
    
        del (essentials_list, feature, gene_alias_dict)
        ##############################################################################
    
        feature_name_list = []
        for feature_name in feature_NameAndType_list:
            feature_name_list.append(feature_name[0])
    
        all_features = {'Feature': feature_name_list,
                        'Essentiality': essentiality_list,
                        'Nreads':N_reads_list,
                        'Ninsertions':N_insrt_list,
                        'Nbasepairs':N_bp_list,
                        'Nreadsperbp':N_reads_per_bp_list,
                        'Ninsertionsperbp':N_insrt_per_bp_list}
    
        if file_counter == 0:
            dna_df2 = pd.DataFrame(all_features, columns = [column_name for column_name in all_features]) #search for feature using: dna_df2.loc[dna_df2['Feature'] == 'CDC42']
        else:
            dna_df2["Ninsertionsperbp_2"] = N_insrt_per_bp_list
    
    
        del (feature_NameAndType_list, feature_name_list, feature_name, f_type, f_previous, f_current, N_reads, N_reads_list, N_insrt, N_insrt_list, N_bp, N_bp_list, bp, i, N_reads_per_bp_list, N_insrt_per_bp_list, all_features)
    
#%% CREATE BAR PLOT TO SHOW DIFFERENCES BETWEEN THE NUMBER OF INSERTIONS PER BP PER REGION
    noncoding_color = '#00918F'#"#003231"
    essential_color = '#00918F'#"#00F28E"
    nonessential_color = '#00918F'#"#F20064"
    codingdna_color = '#00918F'#'#00918f'
#    textcolor = "#003231"
    textsize = 14


    feature_middle_pos_list = []
    sum_bp = 0
    for x in dna_df2['Nbasepairs']:
        feature_middle_pos_list.append(x/2 + sum_bp)
        sum_bp += x
    del (x, sum_bp)
    
    feature_width_list = list(dna_df2['Nbasepairs'])
    

    
    barcolor_list = []
    for feature in dna_df2['Feature']:
        if feature == 'noncoding':
            barcolor_list.append(noncoding_color)
        elif dna_df2.loc[dna_df2['Feature'] == feature]['Essentiality'].iloc[0] == False:
            barcolor_list.append(nonessential_color)
        elif dna_df2.loc[dna_df2['Feature'] == feature]['Essentiality'].iloc[0] == True:
            barcolor_list.append(essential_color)
        elif dna_df2.loc[dna_df2['Feature'] == feature]['Essentiality'].iloc[0] == None:
            barcolor_list.append(codingdna_color)
    del (feature)
    


    ###PLOTTING
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.01)


    #DETERMINE DIFFERENCE BETWEEN THE TWO DATASETS
    Ninsertionsperbp_0_list = list(dna_df2['Ninsertionsperbp'])
    Ninsertionsperbp_1_list = list(dna_df2['Ninsertionsperbp_2'])
    Ninsertionsperbp_d_list = [(Ninsertionsperbp_0_list[i] - Ninsertionsperbp_1_list[i])/max(Ninsertionsperbp_0_list) for i in range(len(Ninsertionsperbp_0_list))]


    ax = plt.subplot(grid[0,0])
    ax.bar(feature_middle_pos_list, Ninsertionsperbp_d_list, feature_width_list, color=barcolor_list)
    if region_start != None and region_end != None and region_start < len_chr and region_end < len_chr:
        ax.set_xlim(region_start, region_end)
    else:
        ax.set_xlim(0, len_chr)
#    ax.set_ylim(0, max(dna_df2['Ninsertionsperbp']) + 0.1*max(dna_df2['Ninsertionsperbp']))
    ax.grid(linestyle='-', alpha=1.0)
    ax.tick_params(labelsize=textsize)
    ax.tick_params(axis='x', which='major')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.xaxis.get_offset_text().set_fontsize(textsize)
    ax.set_xlabel("Basepair position on chromosome "+chrom, fontsize=textsize)
    ax.set_ylabel("Relative difference insertions per bp per region",fontsize=textsize)
    ax.set_title(r"(WT1 - $\Delta$Dpl1) / max(WT1)", fontsize=textsize)
#    ax.set_title(r"(WT2 - WT1 / max(WT2)", fontsize=textsize)
    legend_noncoding = mpatches.Patch(color=noncoding_color, label="Noncoding DNA")
    legend_essential = mpatches.Patch(color=essential_color, label="Annotated essential genes")
    legend_nonessential = mpatches.Patch(color=nonessential_color, label="Non-essential genes")
    legend_coding = mpatches.Patch(color=codingdna_color, label="Other genomic regions")
#    ax.legend(handles=[legend_noncoding, legend_essential, legend_nonessential, legend_coding]) #ADD
    
    
#%% RETURN STATEMENT
    return(dna_df2)
    
#%% Feature position function

def feature_position(feature_dict, chrom, start_chr, dna_dict, feature_type=None):
    
    position_dict = {}
    for feat in feature_dict:
        if feature_dict.get(feat)[5] == chrom:
            if feat.startswith("TEL") and feat.endswith('L'): #correct for the fact that telomeres at the end of a chromosome are stored in the reverse order.
                position_dict[feat] = [feature_dict.get(feat)[5], feature_dict.get(feat)[7], feature_dict.get(feat)[6]]
            else:
                position_dict[feat] = [feature_dict.get(feat)[5], feature_dict.get(feat)[6], feature_dict.get(feat)[7]]


    for feat in position_dict:
        for bp in range(int(position_dict.get(feat)[1])+start_chr, int(position_dict.get(feat)[2])+start_chr):
            if dna_dict[bp] == ['noncoding', None]:
                dna_dict[bp] = [feat, feature_type]
            else:
#                print('Bp %i is already occupied by %s' % (bp, str(dna_dict.get(bp))))
                pass


    return(dna_dict)


#%% executing the function

if __name__ == '__main__':
<<<<<<< HEAD
    dna_df2 = dna_features_compare(region = "III",
                 wig_file_list = [r"N:\tnw\BN\LL\Shared\Gregory\testing_site\Benoit_test_data\wt1_KornmannLab_20200812\ERR1533148_trimmed.sorted.bam.wig",
                             r"N:\tnw\BN\LL\Shared\Gregory\testing_site\Benoit_test_data\dDpl1_KornmannLab_20200803\E-MTAB-4885.Dpl1Kan.sorted.bam.wig"],
                 pergene_insertions_file_list = [r"N:\tnw\BN\LL\Shared\Gregory\testing_site\Benoit_test_data\wt1_KornmannLab_20200812\ERR1533148_trimmed.sorted.bam_pergene_insertions.txt",
                                            r"N:\tnw\BN\LL\Shared\Gregory\testing_site\Benoit_test_data\dDpl1_KornmannLab_20200803\E-MTAB-4885.Dpl1Kan.sorted.bam_pergene_insertions.txt"])
=======
    dna_df2 = dna_features_compare(region = "VIII",
<<<<<<< HEAD
                 wig_file_list = [r"C:\Users\gregoryvanbeek\Documents\testing_site\wt2_testfolder\align_out\ERR1533148_trimmed.sorted.bam.wig",
                             r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig"],
                 pergene_insertions_file_list = [r"C:\Users\gregoryvanbeek\Documents\testing_site\wt2_testfolder\align_out\ERR1533148_trimmed.sorted.bam_pergene_insertions.txt",
                                            r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt"])
>>>>>>> 45dad1c4f93756d74aa008dc8c863c549363fd95
=======
                 wig_file_list = [r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.wig",
                                  r"C:\Users\gregoryvanbeek\Documents\testing_site\dDpl1_testfolder\align_out\E-MTAB-4885.Dpl1Kan.sorted.bam.wig"],
                 pergene_insertions_file_list = [r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam_pergene_insertions.txt",
                                                 r"C:\Users\gregoryvanbeek\Documents\testing_site\dDpl1_testfolder\align_out\E-MTAB-4885.Dpl1Kan.sorted.bam_pergene_insertions.txt"])
>>>>>>> 5183930047024746128f4e376e4177d31be8e315




