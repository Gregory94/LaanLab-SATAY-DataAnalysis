# -*- coding: utf-8 -*-
"""The module includes two functions that have the same purpose, to make a profile plot for a specified chromosome.
The transposon_profile function plots a bar plot for the number of transposons in the chromosome.
The read_profile function plots a bar plot for the number of reads in the chromosome.
The background of the barplots are color coded. A red area indicates a gene that is not annotated as being essential (in a WT background). A green area indicates an annotated essential gene.
Both functions require the modules chromosome_and_gene_positions.py, essential_genes_names.py and gene_names.py including the required files for the functions (see the help in these functions).
"""
#%%
import os, sys
import numpy as np
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic, gene_position
from essential_genes_names import list_known_essentials
from gene_names import gene_aliases
from chromosome_names_in_files import chromosome_name_bedfile, chromosome_name_wigfile

#%%

def transposon_profile(chrom='I',bar_width=None,bed_file = None):
    '''This function creates a bar plot along a specified chromosome for the number of transposons.
    The height of each bar represents the number of transposons at the genomic position indicated on the x-axis.
    The input is as follows: which chromosome (indicated by roman numeral), bar_width, bed_file.
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    The bed_file is one of the files created by the Matlab code from the kornmann-lab.
    The background of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases')
    '''
    #bed_file = r'X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_Trimmed_Aligned\Cerevisiae_WT1_Michel2017_Trimmed_Aligned.sorted.bam.bed'
#%% USED FILES
    gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    essential_genes_files = [os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                            os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_2.txt')]
    gene_information_file = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')
#%% GET CHROMOSOME LENGTHS AND POSITIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    
#%% CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
#    arabic_to_roman_dict, roman_to_arabic_dict = chromosomename_roman_to_arabic()    
#    chromosome_romannames_list = []
#    for roman in roman_to_arabic_dict:
#        chromosome_romannames_list.append(roman)
#    
    
#    for chrom in chromosomenames_list:
#    chrom_index = chromosome_romannames_list.index(chrom)
    print('Chromosome length: ',chr_length_dict.get(chrom))
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/400)
    
    
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]
    genes_essential_list = list_known_essentials(essential_genes_files)
    gene_alias_list = gene_aliases(gene_information_file)[0]
    
#%% READ BED FILE
    with open(bed_file) as f:
        lines = f.readlines()
    
#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
#    chrom_names_dict = {}
#    chrom_start_index_dict = {}
#    chrom_end_index_dict = {}
#    
#    chrom_name = ''
#    chr_counter = 0
#    line_counter = 0
#    stop_loop = False
##    for line in lines:
#    while stop_loop is False:
#        line = lines[line_counter]
#        chrom_name_current = line.split(' ')[0]
#        if not chrom_name_current.startswith('track'): #SKIP HEADER
#            if not chrom_name_current.startswith('chrM'): #SKIP MITOCHRONDRIAL CHROMOSOMES
#                if chrom_name_current != chrom_name:
#                    chrom_names_dict[chromosome_romannames_list[chr_counter]] = chrom_name_current
#                    chrom_name = chrom_name_current
#                    print('Chromosome ',chromosome_romannames_list[chr_counter], 'is ',chrom_name_current)
#                    
#                    chrom_start_index_dict[chromosome_romannames_list[chr_counter]] = line_counter #GET START INDEX IN THE BED FILE OF THE CURENT CHROMOSOME
#                    if chr_counter != 0:
#                        chrom_end_index_dict[chromosome_romannames_list[chr_counter-1]] = line_counter-1 #GET THE END INDEX IN THE BED OF THE PREVIOUS CHROMOSOME (SKIP FOR THE FIRST CHROMOSOME)
#
#                    chr_counter += 1
#
#            elif chrom_name_current.startswith('chrM'):
#                chrom_end_index_dict[chromosome_romannames_list[-1]] = line_counter-1 #GET THE END INDEX IN THE BED FILE FOR THE FINAL CHROMOSOME
#                stop_loop = True
#                
#        line_counter += 1

    chrom_start_index_dict, chrom_end_index_dict= chromosome_name_bedfile(lines)[1:3]

#%% GET ALL TRANSPOSON COUNTS
#    allinsertionsites_list = list(range(0,chr_length_dict.get(chrom)))
    alltransposoncounts_list = np.zeros(chr_length_dict.get(chrom))
    for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
        line = line.strip('\n').split()
        alltransposoncounts_list[int(line[1])] += 1
    
    
    
#%% BINNING OF THE READS
    #THE LIST WITH ALL THE TRANPOSONS FOR THE CURRENT CHROMOSOME IS TYPICALLY REALLY LARGE.
    #TO COMPRESS THIS LIST, THE BASEPAIR POSITIONS ARE GROUPED IN GROUPS WITH SIZE DEFINED BY 'BAR_WIDTH'
    #IN EACH GROUP THE NUMBER OF readS ARE SUMMED UP.
    #THIS IS DONE TO SPEED UP THE SCRIPT AS PLOTTING ALL VALUES IS SLOW
    alltransposoncounts_binnedlist = []
    val_counter = 0
    sum_values = 0
    if bar_width == 1:
        alltransposoncounts_binnedlist = alltransposoncounts_list
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/float(bar_width)))
    else:
        for n in range(len(alltransposoncounts_list)):
            if val_counter % bar_width != 0:
                sum_values += alltransposoncounts_list[n]
            elif val_counter % bar_width == 0:
                alltransposoncounts_binnedlist.append(sum_values)
                sum_values = 0
            val_counter += 1
            
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/bar_width)+1)
    
    
#%% PLOTTING
    print('Plotting chromosome ', chrom, '...')
    print('bar width for plotting is ',bar_width)
    
    textsize = 20
    textcolor = "#000000"

    plt.figure(figsize=(17,6))
    grid = plt.GridSpec(20, 1, wspace=0.0, hspace=0.0)
    
    binsize = bar_width
    ax = plt.subplot(grid[0:19,0])
#    for gene in genes_currentchrom_pos_list:
#        gene_start_pos = int(gene_pos_dict.get(gene)[1])
#        gene_end_pos = int(gene_pos_dict.get(gene)[2])
#        if gene in genes_essential_list:
#            ax.axvspan(gene_start_pos,gene_end_pos,facecolor="#BBE6AA",alpha=1.0)
##            ax.text(gene_start_pos,max(alltransposoncounts_binnedlist),gene_alias_list.get(gene)[0], rotation=90, fontsize=18)
#        else:
#            ax.axvspan(gene_start_pos,gene_end_pos,facecolor="#F6A089",alpha=1.0)
    ax.bar(allinsertionsites_list,alltransposoncounts_binnedlist,width=binsize,color="#00918f")
    ax.tick_params(axis='both', which='major', labelsize=textsize)
    ax.set_axisbelow(True)
    ax.grid(True)
    ax.set_xlim(0,chr_length_dict.get(chrom))
    ax.tick_params(axis='x', which='major', pad=30)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.xaxis.get_offset_text().set_fontsize(textsize)
    ax.set_xlabel("Basepair position on chromosome "+chrom, fontsize=textsize, color=textcolor, labelpad=10)
    ax.set_ylabel('Tranposon count', fontsize=textsize, color=textcolor)
#    ax.set_title('Transposon profile for chromosome '+chrom)


    axc = plt.subplot(grid[19,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            axc.axvspan(gene_start_pos,gene_end_pos,facecolor="#00F28E",alpha=0.8)
#            ax.text(gene_start_pos,max(alltransposoncounts_binnedlist),gene_alias_list.get(gene)[0], rotation=90, fontsize=18)
        else:
            axc.axvspan(gene_start_pos,gene_end_pos,facecolor="#F20064",alpha=0.8)    
    axc.set_xlim(0,chr_length_dict.get(chrom))
    axc.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off

    axc.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left=False,        # ticks along the bottom edge are off
        right=False,       # ticks along the top edge are off
        labelleft=False)   # labels along the bottom edge are off

    plt.show()




#%%
def read_profile(chrom='I',bar_width=None,wig_file = None):
    '''This function creates a bar plot along a specified chromosome for the number of reads.
    The height of each bar represents the number of reads at the genomic position indicated on the x-axis.
    The input is as follows: which chromosome (indicated by roman numeral), bar_width, wig_file.
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    The wig_file is one of the files created by the Matlab code from the kornmann-lab.
    The background of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases')
    '''

#%% USED FILES
    gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    essential_genes_files = [os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                            os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_2.txt')]
    gene_information_file = os.path.join(file_dirname,'Data_Files','Yeast_Protein_Names.txt')
#%%
    #GET CHROMOSOME LENGTHS AND POSITIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    
    
    #CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
    arabic_to_roman_dict, roman_to_arabic_dict = chromosomename_roman_to_arabic()    
    chromosomenames_list = []
    for roman in roman_to_arabic_dict:
        chromosomenames_list.append(roman)
        
#%%
#    chrom_index = chromosomenames_list.index(chrom)
    chrom = chrom.upper()
    print('Chromosome length: ',chr_length_dict.get(chrom))
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/600)
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]
    genes_essential_list = list_known_essentials(essential_genes_files)
    gene_alias_list = gene_aliases(gene_information_file)[0]

#%%
    with open(wig_file) as f:
        lines = f.readlines()

#%% GET THE NAMES OF THE CHROMOSOMES AS USED IN THE WIG FILE.
#    chrom_names_dict = {}
#    chrom_names_counter = 0
#    for line in lines:
#        line.strip('\n')
#        chrom_line = 'variableStep'
#        line_split = line.split(' ')
#        if line_split[0] == chrom_line:
#            chromosome_name_wigfile = line_split[1].replace('chrom=chr','').strip('\n')
#            chrom_names_dict[chromosomenames_list[chrom_names_counter]] = chromosome_name_wigfile
#            print('Chromosome ',chromosomenames_list[chrom_names_counter], 'is ',chromosome_name_wigfile)
#            
#            chrom_names_counter += 1

    chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(lines)

#%% GET ALL LINES WITH THE READS FOR THE CURRENT CHROMOSOME
#    line_counter = 0
#    for line in lines:
#        line = line.strip('\n')
#        if line.endswith('chrom=chr'+chrom_names_dict.get(chromosomenames_list[chrom_index])):
#            wigfile_start_index = line_counter + 1
#        elif chrom_names_dict.get(chrom) == chrom_names_dict.get('XVI'): #CHECK IF THE LAST CHROMOSOME IS REACHED, SINCE THEN THE NEXT CHROMOSOME DOES NOT NEED TO BE SEARCHED AS THIS WON'T EXISTS
#            wigfile_end_index = len(lines)-1 #GET INDEX LAST ELEMENT
#        elif line.endswith('chrom=chr'+chrom_names_dict.get(chromosomenames_list[chrom_index+1])):
#            wigfile_end_index = line_counter
#        line_counter += 1

    wigfile_start_index = chrom_start_line_dict.get(chrom)
    wigfile_end_index = chrom_end_line_dict.get(chrom)

#%%    



#    allinsertionsites_list = list(range(0,chr_length_dict.get(chrom))) #CREATE LIST OF ALL POSIBLE INSERTION SITES IN THE CURRENT CHROMOSOME
    allreadscounts_list = np.zeros(chr_length_dict.get(chrom)) #FOR EACH INSERTION SITE LIST THE NUMBER OF read INSERTION. BY DEFAULT THIS 0 AND IS LATER UPDATED IF AN INSERTION SITE IS PRESENT IN THE WIG FILE
    #GET ALL read COUNTS FOR THE CURRENT CHROMOSOME
    for line in lines[wigfile_start_index:wigfile_end_index]:
        line = line.strip(' \n').split()
        allreadscounts_list[int(line[0])] = int(line[1])
    
#%%    
    
    
    #THE LIST WITH ALL THE TRANPOSONS FOR THE CURRENT CHROMOSOME IS TYPICALLY REALLY LARGE.
    #TO COMPRESS THIS LIST, THE BASEPAIR POSITIONS ARE GROUPED IN GROUPS WITH SIZE DEFINED BY 'BAR_WIDTH'
    #IN EACH GROUP THE NUMBER OF readS ARE SUMMED UP.
    #THIS IS DONE TO SPEED UP THE SCRIPT AS PLOTTING ALL VALUES IS SLOW
#    bar_width = 1000
    allreadscounts_binnedlist = []
    val_counter = 0
    sum_values = 0
    if bar_width == 1:
        allreadscounts_binnedlist = allreadscounts_list
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/float(bar_width)))
    else:
        for n in range(len(allreadscounts_list)):
            if val_counter % bar_width != 0:
                sum_values += allreadscounts_list[n]
            elif val_counter % bar_width == 0:
                allreadscounts_binnedlist.append(sum_values)
                sum_values = 0
            val_counter += 1
            
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/bar_width)+1)

#%%
# =============================================================================
#USE alltransposoncounts_binnedlist FOR NORMALIZING allreadscounts_binnedlist WITH THE AMOUNT OF TRANSPOSONS IN EACH BIN.
#     allreadscounts_binnedlist_np = np.array(allreadscounts_binnedlist, dtype=np.float)
#     alltransposoncounts_binnedlist_np = np.array(alltransposoncounts_binnedlist, dtype=np.float)
#     allreadscounts_binnedlist_np_norm = list(allreadscounts_binnedlist_np/alltransposoncounts_binnedlist_np)
#     allreadscounts_binnedlist_np_norm = np.where(allreadscounts_binnedlist_np_norm==np.inf,0,allreadscounts_binnedlist_np_norm)
# =============================================================================
#%%

    print('Plotting chromosome ', chrom, '...')
    print('bar width for plotting is ',bar_width)

    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(1, 1, wspace=0.0, hspace=0.0)

    textsize = 20

    binsize = bar_width
    ax = plt.subplot(grid[0,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            ax.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
            ax.text(gene_start_pos,max(allreadscounts_binnedlist),gene_alias_list.get(gene)[0], rotation=45)
        else:
            ax.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)
    ax.bar(allinsertionsites_list,allreadscounts_binnedlist,width=binsize,color=[0.0,0.0,0.0,0.8])
    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(True)
    ax.set_xlim(0,chr_length_dict.get(chrom))
    ax.set_xlabel('Basepair position on chromosome '+chrom, fontsize=textsize)
    ax.set_ylabel('Read count (log_10)', fontsize=textsize)
#    ax.set_title('Read profile for chromosome '+chrom)
    plt.tight_layout()



#%%    





#%%
if __name__ == '__main__':
#    read_profile(chrom='xv',wig_file=r"C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\satay_analysis_testdata\Output_Processing\Cerevisiae_WT2_Michel2017_trimmed1.bam.wig")
    transposon_profile(chrom='IX', bar_width=None, bed_file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed")