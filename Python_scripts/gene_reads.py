# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 14:10:52 2020

@author: gregoryvanbeek
"""
#%%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import gene_position
from gene_names import gene_aliases
from chromosome_names_in_files import chromosome_name_bedfile
import statistics_perchromosome

#%%

def gene_reads(gene_name=None,region=None,bed_file=None,savefigure=False):
    '''This script makes a profile plot for the number of reads per tranposon for a specific genomic region.
    Input is a region and the .bed file from the output of the Matlab code from the Kornmann-lab.
    The region can be defined either as a gene name (e.g. 'bem1') or as a list consisting of three elements where the first element is the chromosome name, the start and end position respectively (e.g. ['I',1,4000]).
    If a gene name is input, the script searches in a .gff file (downloaded from yeastgenome.org).
    The output is a bar plot where the number of reads divided by the number of transposons.
    '''
#%% USED FILES
    gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    gene_information_file = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')

#%%SAVE FILES
    if savefigure == True:
        save_figure_path = file_dirname
#%% GET START AND END POSITION OF GENE
    if gene_name.upper() == 'HOLOCUS' or gene_name == 'HO-LOCUS':
        gene_pos = ['IV',46271,48031]
        gene_name = 'HOlocus'

    elif gene_name != None:
        gene_pos_dict = gene_position(gff_file) #GET POSITION INFORMATION OF ALL GENES
        
        gene_name = gene_name.upper() #CAPITALIZE GENE NAME
        if gene_pos_dict.get(gene_name) == None: #CHECK IF GENE_NAME EXISTS IN GENE_POS_DICT. IF NOT, CHECK IF ANY OF THE ALIASES EXISTS
            gene_alias_dict = gene_aliases(gene_information_file)[0]
            gene_alias_key = [k for k, v in gene_alias_dict.items() if gene_name in v]
            print('gene_alias_key ',gene_alias_key[0])
            if gene_pos_dict.get(gene_alias_key[0]) == None: #IF KEY DOES ALSO NOT EXISTS IN GENE_POS_DICT, CHECK IF MORE ALIASES EXISTS OF GENE_NAME
                gene_alias_list = gene_alias_dict.get(gene_alias_key[0])
                for gene_alias in gene_alias_list:
                    if gene_pos_dict.get(gene_alias) != None:
                        gene_pos = gene_pos_dict.get(gene_alias)
                        print('The alias ',gene_alias, ' is used for the gene ',gene_name)
            else:
                gene_pos = gene_pos_dict.get(gene_alias_key[0])
                print('The alias ',gene_alias_key[0], ' is used for the gene ',gene_name)
                
        else:
            gene_pos = gene_pos_dict.get(gene_name)

        gene_orien = gene_pos[3]

    elif region != None:
        gene_pos = region
        

    gene_chr = gene_pos[0]
    gene_start = int(gene_pos[1])
    gene_end = int(gene_pos[2])
    if gene_name != None:
        print(gene_name, ' starts at basepair ',gene_start, ' and ends at basepair ',gene_end, ' in chromosome',gene_chr)
    else:
        print('Selected region starts at basepair ',gene_start, ' and ends at basepair ',gene_end, ' in chromosome',gene_chr)


#%% READ THE BED FILE
    with open(bed_file) as f:
        lines = f.readlines()

#%% GET POSITION FOR THE CHROMOSOMES IN THE BED FILE

    chrom_start_line_dict, chrom_end_line_dict= chromosome_name_bedfile(lines)[1:3]

#%% GET ALL READS WITHIN THE GENE
    insertion_list = []
    read_list = []
    for line in lines[chrom_start_line_dict.get(gene_chr):chrom_end_line_dict.get(gene_chr)]:
        line_list = line.strip('\n').split()
        if gene_start <= int(line_list[1]) <= gene_end:
            insertion_list.append(int(line_list[1]))

            read_value = (int(line_list[4])-100)/20 #the matlab script by benoit takes the number of reads*20+100. This line makes this undone
            read_list.append(read_value)

#%% ACCOUNT FOR DOUBLE INSERTIONS FOR PLOTTING
#see for example chromosome I, bp 3891
    unique_insertion_list = []
    duplicate_insertion_list = []
    for ins in insertion_list: #FIND THE CHROMOSOME POSITION OF ALL DUPLICATED TRANSPOSON INSERTION SITES
        if ins not in unique_insertion_list:
            unique_insertion_list.append(ins)
        else:
            duplicate_insertion_list.append(ins)
    duplicate_insertion_list = np.unique(duplicate_insertion_list) #ACCOUNT FOR THE SITUATION WHERE THERE ARE MORE THAN TWO INSERTIONS AT THE SAME LOCATION

    duplicate_index_list = []
    for dup in duplicate_insertion_list:
        insertion_arr = np.asarray(insertion_list)
        duplicate_index_list.append(np.where(insertion_arr == dup)) #GET ALL INDICES OF THE LIST OF TRANSPOSON INSERTIONS WHERE THE DUPLICATES ARE PRESENT. EACH INSERTION LOCATION GETS ITS OWN NUMPY ARRAY WITHIN THIS LIST

    if len(duplicate_index_list) > 0:
        number_of_duplicates_list = [1]*len(insertion_list) #MAKE LIST OF ONES WITH SAME LENGTH AS INSERTION_LIST FOR STORING NUMBER OF DUPLICATES
        delete_index = []
        for ind_arr in duplicate_index_list: #LOOP OVER ALL INDICES OF DUPLICATES
            ind_list = ind_arr[0]
            ind_list_max = max(ind_list) #GET THE LAST INDEX OF THE DUPLICATES
#            print('Mulitple transposons found at ',ind_list)
            for ind in ind_list:
                if not ind == ind_list_max:
                    read_list[ind_list_max] += read_list[ind] #ADD UP THE READS TO THE LAST DUPLICATE
                    number_of_duplicates_list[ind_list_max] = len(ind_list) #UPDATE NUMBER OF DUPLICATES
                    delete_index.append(ind)

        #REVERSE LOOP OVER LIST FOR DELETING
        for del_ind in reversed(delete_index):
            del read_list[del_ind] #DELETES THE INDEX WHICH IS NOW ADDED UP TO THE LAST INDEX
            del insertion_list[del_ind] #DELETES THE SAME INDICES IN THE INSERTION LISTS.
            del number_of_duplicates_list[del_ind]

        readspertransposon_list = [x/y for x, y in zip(read_list, number_of_duplicates_list)] #DIVIDE THE NUMBER OF READS BY THE NUMBER OF TRANSPOSONS
    else:
        readspertransposon_list = read_list

#%% MAKE LIST OF ALL LOCATIONS IN THE GENE WITH THE NUMBER OF READS IN EACH LOCATION
    gene_length = gene_end-gene_start
    print('Length of region of interest is ',gene_length)
    insertion_roi_list = list(range(gene_start,gene_end+1))
    reads_roi_list = list(np.zeros(gene_length+1))
    
    read_index = 0
    for position in insertion_list:
        roi_index = insertion_roi_list.index(position)
        reads_roi_list[roi_index] = float(readspertransposon_list[read_index])
        read_index += 1

#%% CALCULATE SOME STATISTICAL VALUES FOR THE SELECTED REGION
#insertion_roi_list := list of all potential insertion sites in the region
#reads_roi_list := number of reads in the selected region.

    bp_between_tn_insertions_chr_dict = statistics_perchromosome.chromosome_insertion_periodicity(gene_chr,bed_file)
    insertion_chromosome_avgperiodicity = np.nanmean(bp_between_tn_insertions_chr_dict.get(gene_chr))
    insertion_chromosome_stdperiodicity = np.nanstd(bp_between_tn_insertions_chr_dict.get(gene_chr))
    insertion_chromosome_firstquartileperiodicity = np.nanpercentile(bp_between_tn_insertions_chr_dict.get(gene_chr),25)
    insertion_chromosome_medperiodicity = np.nanpercentile(bp_between_tn_insertions_chr_dict.get(gene_chr),50)
    insertion_chromosome_thirdquartileperiodicity = np.nanpercentile(bp_between_tn_insertions_chr_dict.get(gene_chr),75)

    coverage_percentage = (len(read_list)/len(insertion_roi_list)*100)

    if insertion_list != []:
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(insertion_list[:-1], insertion_list[1:])]
        bp_between_tn_insertions.insert(0,insertion_list[0] - gene_start) #ADD START OF GENE (bp=0)
        bp_between_tn_insertions.append(gene_end - insertion_list[-1]) #ADD END OF GENE (bp=INDEX LAST TN - GENE LENGTH)

        max_empty_region = max(bp_between_tn_insertions)
        
        insertion_avgperiodicity = np.nanmean(bp_between_tn_insertions)
        insertion_stdperiodicity = np.nanstd(bp_between_tn_insertions)
        insertion_firstquartileperiodicity = np.nanpercentile(bp_between_tn_insertions,25)
        insertion_medperiodicity = np.nanpercentile(bp_between_tn_insertions,50)
        insertion_thirdquartileperiodicity = np.nanpercentile(bp_between_tn_insertions,75)
    else:
        max_empty_region = gene_length
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(insertion_list[:-1], insertion_list[1:])]
        insertion_avgperiodicity = 0
        insertion_medperiodicity = 0

    if insertion_list != []:
        print('')
        print('Percentage of coverage is %.2f' % coverage_percentage)
        print('Largest area without insertions is %.0f' % max_empty_region)
        print('')
        print('Mean transposon insertion frequency in gene is %.2f, %.2f ' % (insertion_avgperiodicity, insertion_stdperiodicity))
        print('Mean transposon insertion frequency in chromosome is %.2f, %.2f' % (insertion_chromosome_avgperiodicity, insertion_chromosome_stdperiodicity))
        print('')
        print('Quartiles transposon insertion frequency in gene is %.2f, %.2f, %.2f ' % (insertion_firstquartileperiodicity, insertion_medperiodicity, insertion_thirdquartileperiodicity))
        print('Quartiles transposon insertion frequency in chromosome is %.2f, %.2f, %.2f ' % (insertion_chromosome_firstquartileperiodicity, insertion_chromosome_medperiodicity, insertion_chromosome_thirdquartileperiodicity))


#%% BINNING OF THE READS_ROI_LIST
### DETERMINE VARIABLE BIN WIDTH SUCH THAT EITHER EACH BIN INCLUDES 8 TN OR THE WIDTH OF THE BAR IS EIGTH TIMES THE AVERAGE DISTANCE BETWEEN TN IN THE CHROMOSOME.

#    bin_width = int(insertion_chromosome_avgperiodicity*8)
#    reads_roi_binnedlist = []
#    val_counter = 0
#    sum_values = 0
#
#    for n in range(len(reads_roi_list)):
#        if val_counter % bin_width != 0:
#            sum_values += reads_roi_list[n]
#        elif val_counter % bin_width == 0:
#            reads_roi_binnedlist.append(sum_values)
#            sum_values = 0
#        val_counter += 1
#
#    insertion_roi_binnedlist = np.linspace(gene_start,gene_end,int(gene_length/bin_width)+1)

    reads_roi_binnedlist = [] #STORES THE HEIGHT OF THE BARS
    insertion_roi_binnedlist = [] #STORES THE X POSITION OF THE BARS
    bin_width = [] #STORES THE WIDTH OF THE BARS

    currentbin_reads_list = []
    currentbin_insertion_list = []
    reads_currentbin = 0
    insertion_counter = 0
    for ins in range(0,len(insertion_roi_list)):
        currentbin_reads_list.append(reads_roi_list[ins])
        currentbin_insertion_list.append(insertion_roi_list[ins])
        if reads_roi_list[ins] > 0:
            reads_currentbin += 1

        if reads_currentbin == 8: #STOP IF THERE ARE 8 INSERTIONS FOUND
            reads_roi_binnedlist.append(sum(currentbin_reads_list))
            bin_width.append(max(currentbin_insertion_list) - min(currentbin_insertion_list))
            insertion_roi_binnedlist.append(insertion_roi_list[ins]-bin_width[-1]/2)
            
            currentbin_reads_list = []
            currentbin_insertion_list = []
            reads_currentbin = 0
            insertion_counter = 0

        elif insertion_counter == int(insertion_chromosome_avgperiodicity*8): #STOP IF THE LENGTH OF THE CURRENTBIN EXCEEDS AVERAGE PERIODICITY OF THE CHROMOSOME * 8
            reads_roi_binnedlist.append(sum(currentbin_reads_list))
            bin_width.append(max(currentbin_insertion_list) - min(currentbin_insertion_list))
            insertion_roi_binnedlist.append(insertion_roi_list[ins]-bin_width[-1]/2)
            
            currentbin_reads_list = []
            currentbin_insertion_list = []
            reads_currentbin = 0
            insertion_counter = 0

        elif gene_start+ins == gene_end: #ACCOUNT FOR THE FINAL BIN IN THE GENE.
            reads_roi_binnedlist.append(sum(currentbin_reads_list))
            bin_width.append(max(currentbin_insertion_list) - min(currentbin_insertion_list))
            insertion_roi_binnedlist.append(insertion_roi_list[ins]-bin_width[-1]/2)

        insertion_counter += 1

#%% MAKE BAR PLOT FOR READS IN CHROMOSOME

    if gene_name != None:
        print('Plotting reads for gene ', gene_name, '...')
    else:
        print('Plotting reads in range ', gene_start, '..', gene_end, 'in chromosome ', gene_chr, '...')

    text_size = 12

    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(1, 3, wspace=0.4, hspace=0.3)
    
    ax = plt.subplot(grid[0,:2])
    ax.bar(insertion_roi_binnedlist,reads_roi_binnedlist,width=bin_width,facecolor=np.array([126.0,164.0,179.0])/255,edgecolor='w')
    ax.set_axisbelow(True)
    ax.grid(True)
    if gene_name != None:
        ax.set_title(gene_name, fontweight='bold', fontsize=text_size)
    elif region == ['IV',46271,48031]:
        ax.set_title('HO-locus', fontweight='bold', fontsize=text_size)
    else:
        ax.set_title(str(gene_chr) + str(gene_start) + '-' + str(gene_end))
    ax.set_xlabel('Basepair position in chromosome '+ gene_chr, fontsize=text_size) 
    ax.set_ylabel('Read/Tn', fontsize=text_size)
    ax.set_xlim(gene_start,gene_end)

    if gene_name != None and gene_name != 'HOlocus':
        if insertion_list != []:
            textstr = '\n'.join((r'Reading orientation of gene: ' + gene_orien,
                                r'Transposon coverage = %.2f percent' % (coverage_percentage),
                                r'Mean transposon insertion frequency in gene is %.2f, %.2f ' % (insertion_avgperiodicity, insertion_stdperiodicity),
                                r'Mean transposon insertion frequency in chromosome is %.2f, %.2f' % (insertion_chromosome_avgperiodicity, insertion_chromosome_stdperiodicity),
                                r'Quartiles transposon insertion frequency in gene is %.2f, %.2f, %.2f ' % (insertion_firstquartileperiodicity, insertion_medperiodicity, insertion_thirdquartileperiodicity),
                                r'Quartiles transposon insertion frequency in chromosome is %.2f, %.2f, %.2f ' % (insertion_chromosome_firstquartileperiodicity, insertion_chromosome_medperiodicity, insertion_chromosome_thirdquartileperiodicity)
                                ))
        else:
            textstr = (r'Reading orientation of gene: ' + gene_orien)
    else:
        textstr = '\n'.join((r'Transposon coverage = %.2f percent' % (coverage_percentage),
                            r'Mean transposon insertion frequency in gene is %.2f, %.2f ' % (insertion_avgperiodicity, insertion_stdperiodicity),
                            r'Mean transposon insertion frequency in chromosome is %.2f, %.2f' % (insertion_chromosome_avgperiodicity, insertion_chromosome_stdperiodicity),
                            r'Quartiles transposon insertion frequency in gene is %.2f, %.2f, %.2f ' % (insertion_firstquartileperiodicity, insertion_medperiodicity, insertion_thirdquartileperiodicity),
                            r'Quartiles transposon insertion frequency in chromosome is %.2f, %.2f, %.2f ' % (insertion_chromosome_firstquartileperiodicity, insertion_chromosome_medperiodicity, insertion_chromosome_thirdquartileperiodicity)
                            ))
    props = dict(boxstyle='round', facecolor='grey', alpha=0.8)
    ax.text(0.05,0.9, textstr, transform=ax.transAxes, fontsize=text_size,
        verticalalignment='top',horizontalalignment='left', bbox=props)

    for ins in insertion_list: #PLOT TICKS ON THE X AXIS TO INDICATE THE EXACT INSERTION LOCATIONS
        ax.axvline(x=ins, ymax=0.05, linewidth=1, color='k')



# COMPARE DISTRIBUTION OF BASEPAIRS BETWEEN INSERTIONS FOR THE CHROMOSOME AND THE GENE USING VIOLINPLOT
    plt.subplot(grid[0,2])
    
    if gene_name == None:
        if region == ['IV',46271,48031]:
            gene_name = 'HOlocus'
        else:
            gene_name = str(gene_chr) + str(gene_start) + '-' + str(gene_end)
    
    bp_between_tn_insertions_dict = {}
    bp_between_tn_insertions_dict[gene_chr] = bp_between_tn_insertions
    
    df_chr = pd.DataFrame(bp_between_tn_insertions_chr_dict)
    df = pd.DataFrame(bp_between_tn_insertions_dict)
    if len(df) !=0:
        df_concat = pd.concat([df,df_chr], axis=0, ignore_index=True)
        names_list = ['gene']*len(bp_between_tn_insertions) + ['chromosome']*len(bp_between_tn_insertions_chr_dict.get(gene_chr))
        df_concat['label'] = names_list
        xlabel=gene_name+' | '+gene_chr
        df_concat[xlabel] = ''
        df_concat.columns = ['bp between tn','label',xlabel]

        sb.set(style="whitegrid", palette="pastel", color_codes=True)
        sb.violinplot(data=df_concat,x=xlabel, y='bp between tn', hue='label', inner='quartile', scale='width', gridsize=5000, split=True, cut=0, palette={"gene": "g", "chromosome": "r"})
        plt.ylim(0,max_empty_region+10)
    else:
        sb.violinplot(data=df_chr, inner='quartile', orien='v', scale='width', gridsize=5000, cut=0)

    plt.show()
    if savefigure == True:
        plt.savefig(os.path.join(save_figure_path, gene_name+'_TnInsertions.png'), dpi=300)
        print('Figure saved at ', save_figure_path)

#%%
if __name__ == '__main__':
#    gene_reads(region=['I',1,4000],bed_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed")
#    gene_reads(region=['IV',46271,48031],bed_file=r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed", savefigure=True)
    gene_reads(gene_name='bem2',bed_file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed", savefigure=False)
