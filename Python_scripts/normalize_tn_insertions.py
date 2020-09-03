# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:02:44 2020

@author: gregoryvanbeek

This scripts takes a user defined genomic region (i.e. chromosome number, region or gene) and determines the number of transposon insertions in the genomic features (i.e. genes, nc-DNA etc.).
This can be used to determine the number of transposon insertions outside the genes to use this for normalization of the number of transposons in the genes.
"""

#%%
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic
from chromosome_names_in_files import chromosome_name_wigfile
from gene_names import list_gene_names, gene_aliases


#%% USER INPUT AND FILES

region = 'III'

wig_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533148_trimmed.sorted.bam.wig"

pergene_insertions_file = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533148_trimmed.sorted.bam_pergene_insertions.txt"

essentials_file = r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\Python_scripts\Data_Files\Cerevisiae_AllEssentialGenes_List.txt"

gene_information_file = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')

gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')


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

#%% READ WIG FILE FOR GETTING ALL TN INSERTIONS

with open(wig_file, 'r') as f:
    lines = f.readlines()

chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(lines)[1:]

insrt_in_chrom_list = []
reads_in_chrom_list = []
for l in lines[chrom_start_line_dict.get(chrom):chrom_end_line_dict.get(chrom)]:
    insrt_in_chrom_list.append(int(l.strip('\n').split(' ')[0]))
    reads_in_chrom_list.append(int(l.strip('\n').split(' ')[1]))


del (lines, l, f, chrom_start_line_dict, chrom_end_line_dict)

#%% READ PERGENE_INSERTIONS FILE FOR LOCATION AND INSERTIONS FOR EACH GENE.

with open(pergene_insertions_file) as f:
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
    
        gene_position_dict[genename] = [gene_chrom, gene_start, gene_end]
    
    
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
for bp in range(start_chr, end_chr + 1):
    dna_dict[bp] = 'noncoding'


for gene in gene_position_dict: #IGNORE GENES THAT END WITH -A AS THESE ARE DUBIOUS (?) 
    for bp in range(gene_position_dict.get(gene)[1]+start_chr, gene_position_dict.get(gene)[2]+start_chr+1):
        dna_dict[bp] = gene
        
# GET FEATURES FROM INTERGENIC REGIONS (-> SEE SGD_features.tab IN DATA_FILES IN GITHUB FOLDER)

del (gene, bp)

#%% DETERMINE THE NUMBER OF TRANSPOSONS PER BP FOR EACH FEATURE AND STORE RESULTS IN DATAFRAME

dna_df = pd.DataFrame(list(dna_dict.items()), columns=["BP", "Feature"])


reads_loc_list = [0] * len(dna_dict) # CONTAINS ALL READS JUST LIKE READS_IN_CHROM_LIST, BUT THIS LIST HAS THE SAME LENGTH AS THE NUMBER OF BP IN THE CHROMOSOME WHERE THE LOCATIONS WITH NO READS ARE FILLED WITH ZEROS
i = 0
for ins in insrt_in_chrom_list:
    reads_loc_list[ins] = reads_in_chrom_list[i]
    i += 1
dna_df["Reads"] = reads_loc_list


del (i, ins)

#%% CREATE DATAFRAME FOR EACH FEATURE (E.G. NONCODING DNA, GENE, ETC.) IN THE CHROMOSOME AND DETERMINE THE NUMBER OF INSERTIONS AND READS PER FEATURE.
feature_name_list = []
f_previous = dna_dict.get(start_chr)
N_reads = 0
N_reads_list = []
N_insrt = 0
N_insrt_list = []
N_bp = 0
N_bp_list = []
i = 0
for bp in dna_dict:
    f_current = dna_dict.get(bp)
    if f_current == f_previous:
        N_bp += 1
        N_reads += reads_loc_list[i]
        if not reads_loc_list[i] == 0:
            N_insrt += 1
    elif (f_current != f_previous or (i+start_chr) == end_chr):# and not f_current.endswith('-A'):
        feature_name_list.append(f_previous)
        N_reads_list.append(N_reads)
        N_insrt_list.append(N_insrt)
        N_bp_list.append(N_bp)
        N_reads = 0
        N_insrt = 0
        N_bp = 0
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

essentiality_list = []
for feature in feature_name_list:
    if not feature == 'noncoding':
        if feature in essentials_list:
            essentiality_list.append(True)
        else:
            essentiality_list.append(False)
    else:
        essentiality_list.append(False)

del (essentials_list, feature, gene_alias_dict)
##############################################################################



all_features = {'Feature': feature_name_list,
                'Essentiality': essentiality_list,
                'Nreads':N_reads_list,
                'Ninsertions':N_insrt_list,
                'Nbasepairs':N_bp_list,
                'Nreadspertn':N_reads_per_bp_list,
                'Ninsertionspertn':N_insrt_per_bp_list}

dna_df2 = pd.DataFrame(all_features, columns = [column_name for column_name in all_features]) #search for feature using: dna_df2.loc[dna_df2['Feature'] == 'CDC42']
#CREATE NEW COLUMN WITH ALL DOMAINS OF THE GENE (IF PRESENT) AND ANOTHER COLUMN THAT INCLUDES LISTS OF THE BP POSITIONS OF THESE DOMAINS



del (feature_name_list, f_previous, f_current, N_reads, N_reads_list, N_insrt, N_insrt_list, N_bp, N_bp_list, bp, i, N_reads_per_bp_list, N_insrt_per_bp_list, all_features)


#%% CREATE BAR PLOT
noncoding_color = "#645f9f"
essential_color = "#70ff96"
nonessential_color = "#FF6161" 
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
del (feature)




###PLOTTING
fig = plt.figure(figsize=(19,9))
grid = plt.GridSpec(20, 1, wspace=0.0, hspace=0.01)


ax = plt.subplot(grid[0:19,0])
ax.bar(feature_middle_pos_list, list(dna_df2['Ninsertionspertn']), feature_width_list, color=barcolor_list)
if region_start != None and region_end != None and region_start < len_chr and region_end < len_chr:
    ax.set_xlim(region_start, region_end)
else:
    ax.set_xlim(0, len_chr)
ax.set_ylim(0, max(dna_df2['Ninsertionspertn']) + 0.1*max(dna_df2['Ninsertionspertn']))
ax.grid(linestyle='--', alpha=0.5)
ax.tick_params(labelsize=textsize)
ax.set_xticklabels([])
legend_noncoding = mpatches.Patch(color=noncoding_color, label="Noncoding DNA")
legend_essential = mpatches.Patch(color=essential_color, label="Annotated essential genes")
legend_nonessential = mpatches.Patch(color=nonessential_color, label="Non-essential genes")
ax.legend(handles=[legend_noncoding, legend_essential, legend_nonessential])


axc = plt.subplot(grid[19,0])
ess_start_pos_list = []
ness_start_pos_list = []
ess_end_pos_list = []
ness_end_pos_list = []
l = 0
counter = 0
for width in feature_width_list:
    if dna_df2.loc[counter][1] == True:
#        ess_start_pos_list.append(l)
#        ess_end_pos_list.append(l - 1)
        axc.axvspan(l,l+width,facecolor=essential_color,alpha=0.3)
    elif dna_df2.loc[counter][1] == False and not dna_df2.loc[counter][0] == 'noncoding':
        axc.axvspan(l,l+width,facecolor=nonessential_color,alpha=0.3)
    l += width
    counter += 1
if region_start != None and region_end != None and region_start < len_chr and region_end < len_chr:
    axc.set_xlim(region_start, region_end)
else:
    axc.set_xlim(0, len_chr)
axc.tick_params(labelsize=textsize)
axc.set_yticklabels([])


#FOLLOWING IS FOR SHOWING COMMON AXIS LABELS.
# add a big axis, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.xlabel("bp position in chromosome " + chrom, fontsize=textsize)
plt.ylabel("Transposons/bp per region", fontsize=textsize)



del (ess_start_pos_list, ness_start_pos_list, ess_end_pos_list, ness_end_pos_list, l, counter, width)

#%%

















