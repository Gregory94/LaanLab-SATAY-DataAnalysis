# -*- coding: utf-8 -*-
"""The module includes two functions that have the same purpose, to make a profile plot for a specified chromosome.
The transposon_profile function plots a bar plot for the number of transposons in the chromosome.
The read_profile function plots a bar plot for the number of reads in the chromosome.
The background of the barplots are color coded. A red area indicates a gene that is not annotated as being essential (in a WT background). A green area indicates an annotated essential gene.
Both functions require the modules chromosome_and_gene_positions.py, essential_genes_names.py and gene_names.py including the required files for the functions (see the help in these functions).
"""
#%%
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(1,r'C:\Users\gregoryvanbeek\Documents\GitHub\LaanLab-SATAY-DataAnalysis\python_modules')
from chromosome_and_gene_positions import chromosome_position, gene_position
from essential_genes_names import list_known_essentials
from gene_names import gene_aliases
from chromosome_names_in_files import chromosome_name_wigfile, chromosome_name_bedfile

#%%

def transposon_profile(chrom='I',bar_width=None,bed_files = None):
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
    gff_file = r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    essential_genes_files = [r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_1.txt',
                            r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_2.txt']
    gene_information_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt'
#%% GET CHROMOSOME LENGTHS AND POSITIONS
    chrom=chrom.upper()

    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    print('Chromosome length: ',chr_length_dict.get(chrom))
    
#%% CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/500)
    
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]
    genes_essential_list = list_known_essentials(essential_genes_files)
    gene_alias_list = gene_aliases(gene_information_file)[0]
    
#%% READ BED FILE
    allinsertionsites_allfiles_list = []
    alltransposoncounts_allfiles_binnedlist = []
    for bed_file in bed_files:
        with open(bed_file) as f:
            lines = f.readlines()
    
#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
        chrom_start_index_dict, chrom_end_index_dict= chromosome_name_bedfile(lines)[1:3]

#%% GET ALL TRANSPOSON COUNTS
        alltransposoncounts_list = np.zeros(chr_length_dict.get(chrom))
        for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
            line = line.strip('\n').split()
            alltransposoncounts_list[int(line[1])] += 1

#%% BINNING OF THE READS
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
    
        allinsertionsites_allfiles_list.append(allinsertionsites_list)
        alltransposoncounts_allfiles_binnedlist.append(alltransposoncounts_binnedlist)

#%% DETERMINE DIFFERENCE BETWEEN DATASETS TRANSPOSONCOUNTS
    transposoncounts_positivedifference_list = [0]*len(alltransposoncounts_allfiles_binnedlist[0])
    transposoncounts_negativedifference_list = [0]*len(alltransposoncounts_allfiles_binnedlist[0])
    for i in range(0,len(alltransposoncounts_allfiles_binnedlist[0])):
        difference = alltransposoncounts_allfiles_binnedlist[0][i]-alltransposoncounts_allfiles_binnedlist[1][i]
        if difference >= 0:
            transposoncounts_positivedifference_list[i] = difference
        elif difference < 0:
            transposoncounts_negativedifference_list[i] = -difference

#%% PLOTTING
    print('Plotting chromosome ', chrom, '...')
    print('bar width for plotting is ',bar_width)
    binsize = bar_width
    max_ylim = max([item for sublist in alltransposoncounts_allfiles_binnedlist for item in sublist]) #GET MAXIMUM VALUE FOR SETTING THE Y AXIS LIMIT EQUAL FOR BOTH GRAPHS
    max_ylim = max_ylim + 0.1*max_ylim
    
    
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)


    ax1 = plt.subplot(grid[0,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
            ax1.text(gene_start_pos,max_ylim,gene_alias_list.get(gene)[0], rotation=45)
        else:
            ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)

    ax1.bar(allinsertionsites_allfiles_list[0],alltransposoncounts_allfiles_binnedlist[0],width=binsize,color=(0.2,0.2,0.2,0.8))
    ax1.bar(allinsertionsites_allfiles_list[0],transposoncounts_positivedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8))

    ax1.set_axisbelow(True)
    ax1.grid(True)
    ax1.set_ylabel('Tranposon count Dataset 1', fontsize=12)
    ax1.set_ylim(0,max_ylim)


    ax2 = plt.subplot(grid[1,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
        else:
            ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)

    ax2.bar(allinsertionsites_allfiles_list[1],alltransposoncounts_allfiles_binnedlist[1],width=binsize,color=(0.2,0.2,0.2,0.8), label='Number of transposons')
    ax2.bar(allinsertionsites_allfiles_list[1],transposoncounts_negativedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8), label='Absolute difference datasets (set1-set2)')

    ax2.set_axisbelow(True)
    ax2.grid(True)
    ax2.set_ylabel('Transposon count Dataset 2', fontsize=12)
    ax2.set_xlabel('Basepair position on chromosome '+chrom, fontsize=12)
    ax2.set_ylim(0,max_ylim)
    ax2.invert_yaxis()
    ax2.legend(loc='lower right')
    
    plt.tight_layout()









#%%
def read_profile(chrom='I',bar_width=None,wig_files = None):
    '''This function creates a bar plot along a specified chromosome for the number of reads.
    The height of each bar represents the number of reads at the genomic position indicated on the x-axis.
    The input is as follows: which chromosome (indicated by roman numeral), bar_width, wig_file.
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    The wig_file is one of the files created by the Matlab code from the kornmann-lab.
    The background of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases')
    '''

#%% USED FILES
    gff_file = r"X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3"
    essential_genes_files = [r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_1.txt',
                            r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Cervisiae_EssentialGenes_List_2.txt']
    gene_information_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt'

#%% GET CHROMOSOME LENGTHS AND POSITIONS
    chrom=chrom.upper()
    
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    print('Chromosome length: ',chr_length_dict.get(chrom))
    
#%% CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/500)
        
#%%
#    chrom_index = chromosomenames_list.index(chrom)
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/500)
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]
    genes_essential_list = list_known_essentials(essential_genes_files)
    gene_alias_list = gene_aliases(gene_information_file)[0]

#%% READ WIG FILE
    allinsertionsites_allfiles_list = []
    allreadscounts_allfiles_binnedlist = []
    for wig_file in wig_files:
        with open(wig_file) as f:
            lines = f.readlines()

#%% GET THE NAMES OF THE CHROMOSOMES AS USED IN THE WIG FILE.
        chrom_names_dict, chrom_start_line_dict, chrom_end_line_dict = chromosome_name_wigfile(lines)

#%% GET ALL LINES WITH THE READS FOR THE CURRENT CHROMOSOME
        wigfile_start_index = chrom_start_line_dict.get(chrom)
        wigfile_end_index = chrom_end_line_dict.get(chrom)

#%% DETERMINE THE NUMBER OF READS OF ALL POSSIBLE INSERTION SITES IN THE CHROMOSOME
        allreadscounts_list = np.zeros(chr_length_dict.get(chrom)) #FOR EACH INSERTION SITE LIST THE NUMBER OF read INSERTION. BY DEFAULT THIS 0 AND IS LATER UPDATED IF AN INSERTION SITE IS PRESENT IN THE WIG FILE
        for line in lines[wigfile_start_index:wigfile_end_index]:
            line = line.strip(' \n').split()
            allreadscounts_list[int(line[0])] = int(line[1])
    
#%% BINNING OF THE INSERTION SITES FOR TO SPEED UP THE PLOTTING PROCESS
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

        allinsertionsites_allfiles_list.append(allinsertionsites_list)
        allreadscounts_allfiles_binnedlist.append(allreadscounts_binnedlist)

#%% DETERMINE DIFFERENCE BETWEEN DATASETS TRANSPOSONCOUNTS
    readcounts_positivedifference_list = [0]*len(allreadscounts_allfiles_binnedlist[0])
    readcounts_negativedifference_list = [0]*len(allreadscounts_allfiles_binnedlist[0])
    for i in range(0,len(allreadscounts_allfiles_binnedlist[0])):
        difference = allreadscounts_allfiles_binnedlist[0][i]-allreadscounts_allfiles_binnedlist[1][i]
        if difference >= 0:
            readcounts_positivedifference_list[i] = difference
        elif difference < 0:
            readcounts_negativedifference_list[i] = -difference

#%% PLOTTING
    print('Plotting chromosome ', chrom, '...')
    print('bar width for plotting is ',bar_width)
    binsize = bar_width
    max_ylim = max([item for sublist in allreadscounts_allfiles_binnedlist for item in sublist]) #GET MAXIMUM VALUE FOR SETTING THE Y AXIS LIMIT EQUAL FOR BOTH GRAPHS
    max_ylim = max_ylim + 0.1*max_ylim
    
    plt.figure(figsize=(19,9))
    grid = plt.GridSpec(2, 1, wspace=0.0, hspace=0.0)


    ax1 = plt.subplot(grid[0,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
            ax1.text(gene_start_pos,max_ylim,gene_alias_list.get(gene)[0], rotation=45)
        else:
            ax1.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)
    
    ax1.bar(allinsertionsites_allfiles_list[0],allreadscounts_allfiles_binnedlist[0],width=binsize,color=[0.0,0.0,0.0,0.6])
    ax1.bar(allinsertionsites_allfiles_list[0],readcounts_positivedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8))

    ax1.set_yscale('log')
    ax1.set_axisbelow(True)
    ax1.grid(True)
    ax1.set_ylabel('Read count (log_10) Dataset 1', fontsize=12)
    ax1.set_xticklabels([])
    ax1.set_ylim(0,max_ylim)


    ax2 = plt.subplot(grid[1,0])
    for gene in genes_currentchrom_pos_list:
        gene_start_pos = int(gene_pos_dict.get(gene)[1])
        gene_end_pos = int(gene_pos_dict.get(gene)[2])
        if gene in genes_essential_list:
            ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='g',alpha=0.3)
        else:
            ax2.axvspan(gene_start_pos,gene_end_pos,facecolor='r',alpha=0.3)
    
    ax2.bar(allinsertionsites_allfiles_list[1],allreadscounts_allfiles_binnedlist[1],width=binsize,color=[0.0,0.0,0.0,0.6], label='Number of transposons')
    ax2.bar(allinsertionsites_allfiles_list[1],readcounts_negativedifference_list,width=binsize,color=(0.52,0.71,0.90,0.8), label='Absolute difference datasets (set1-set2)')

    ax2.set_yscale('log')
    ax2.set_axisbelow(True)
    ax2.grid(True)
    ax2.set_ylabel('Read count (log_10) Dataset 2', fontsize=12)
    ax2.set_xlabel('Basepair position on chromosome '+chrom, fontsize=12)
    ax2.set_ylim(0,max_ylim)
    ax2.invert_yaxis()
    ax2.legend(loc='lower right')

    plt.tight_layout()





#%%
if __name__ == '__main__':
#    read_profile(chrom='i',wig_files=[r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.wig",
#                                      r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_Trimmed_Aligned\Cerevisiae_WT1_Michel2017_Trimmed_Aligned.sorted.bam.wig"])
    transposon_profile(chrom='i',bed_files=[r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_ProcessedByBenoit\E-MTAB-4885.WT1.bam.bed",
                                      r"X:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\Michel2017_WT1_SeqData\Cerevisiae_WT1_Michel2017_Trimmed_Aligned\Cerevisiae_WT1_Michel2017_Trimmed_Aligned.sorted.bam.bed"])