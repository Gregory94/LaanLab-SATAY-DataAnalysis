# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 10:30:46 2021

@author: gregoryvanbeek
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from chromosome_and_gene_positions import chromosome_position, gene_position
from essential_genes_names import list_known_essentials
from chromosome_names_in_files import chromosome_name_bedfile



#%% INPUT

bed_file=r"C:\Users\gregoryvanbeek\Documents\Data_Sets\dataset_leila\dataset_leila_wt\leila_dataset_wt_processing\WT_merged-techrep-a_techrep-b_processing\WT_merged-techrep-a_techrep-b_trimmed.sorted.bam.bed"
variable="reads" #"reads" "transposons"
chrom="III"
bar_width=100
savefig=True

#%%
def profile_plot(bed_file, variable="transposons", chrom='I', bar_width=None, savefig=False):
    '''This function creates a bar plot along a specified chromosome for the number of transposons or reads.
    The height of each bar represents the number of transposons or reads at the genomic position indicated on the x-axis.
    The input is as follows: 
        - bed_file: input absolute path to bed file
        - variable: either transposons or reads
        - chrom: roman numeral indicated the chromosome that needs to be plotted
        - bar_width: integer. By default, the bar_width is set to length_chromosome/800
        - savefig: whether to save the figure at the location of the bed file (True or False)
        
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    The bottom part of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases')
    '''
#%% USED FILES
    gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    essential_genes_files = [os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                            os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_2.txt')]

#%% GET CHROMOSOME LENGTHS AND POSITIONS
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    
#%% CREATE LIST OF ALL CHROMOSOMES IN ROMAN NUMERALS
    print('Chromosome length: ',chr_length_dict.get(chrom))
    if bar_width == None:
        bar_width = int(chr_length_dict.get(chrom)/800)
    
    
#%% GET ALL GENES IN CURRENT CHROMOSOME
    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items() if chrom in v]
    genes_essential_list = list_known_essentials(essential_genes_files)

    
#%% READ BED FILE
    with open(bed_file) as f:
        lines = f.readlines()
    
#%% GET NAMES FOR THE CHROMOSOMES IN THE BED FILE
    chrom_start_index_dict, chrom_end_index_dict= chromosome_name_bedfile(lines)[1:3]

#%% GET ALL TRANSPOSON COUNTS
    allcounts_list = np.zeros(chr_length_dict.get(chrom) + 1)
    if variable == "transposons":
        for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
            line = line.strip('\n').split()
            allcounts_list[int(line[1]) - 1] += 1

    elif variable == "reads":
        for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
            line = line.strip('\n').split()
            allcounts_list[int(line[1]) - 1] += (int(line[4])-100)/20

    else:
        print("ERROR: No valid variable argument given. Use transposons or reads")
        exit(1)
                


#%% BINNING OF THE READS
    #THE LIST WITH ALL THE TRANPOSONS FOR THE CURRENT CHROMOSOME IS TYPICALLY REALLY LARGE.
    #TO COMPRESS THIS LIST, THE BASEPAIR POSITIONS ARE GROUPED IN GROUPS WITH SIZE DEFINED BY 'BAR_WIDTH'
    #IN EACH GROUP THE NUMBER OF readS ARE SUMMED UP.
    #THIS IS DONE TO SPEED UP THE SCRIPT AS PLOTTING ALL VALUES IS SLOW
    allcounts_binnedlist = []
    val_counter = 0
    sum_values = 0
    if bar_width == 1:
        allcounts_binnedlist = allcounts_list
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/float(bar_width)))
    else:
        for n in range(len(allcounts_list)):
            if val_counter % bar_width != 0:
                sum_values += allcounts_list[n]
            elif val_counter % bar_width == 0:
                allcounts_binnedlist.append(sum_values)
                sum_values = 0
            val_counter += 1
            
        allinsertionsites_list = np.linspace(0,chr_length_dict.get(chrom),int(chr_length_dict.get(chrom)/bar_width)+1)
    
    
#%% PLOTTING
    print('Plotting chromosome ', chrom, '...')
    print('bar width for plotting is ',bar_width)
    
    textsize = 18
    textcolor = "#000000"

    plt.figure(figsize=(19,9))#(17,6))
    grid = plt.GridSpec(20, 1, wspace=0.0, hspace=0.0)
    
    binsize = bar_width
    ax = plt.subplot(grid[0:19,0])
    ax.bar(allinsertionsites_list,allcounts_binnedlist,width=binsize,color="#000000")
    ax.tick_params(axis='both', which='major', labelsize=textsize)
    ax.set_axisbelow(True)
    ax.grid(True)
    ax.set_xlim(0,chr_length_dict.get(chrom))
#    ax.set_ylim(0, 200)
    ax.tick_params(axis='x', which='major', pad=30)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.xaxis.get_offset_text().set_fontsize(textsize)
    ax.set_xlabel("Basepair position on chromosome "+chrom, fontsize=textsize, color=textcolor, labelpad=10)
    if variable == "transposons":
        ax.set_ylabel('Transposon count', fontsize=textsize, color=textcolor, labelpad=25)
    elif variable == "reads":
        ax.set_ylabel('Read count', fontsize=textsize, color=textcolor, labelpad=25)
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

    if savefig == True and variable == "transposons":
        savepath = os.path.splitext(bed_file)
        print('saving figure at %s' % savepath[0]+'_transposonplot_chrom'+chrom+'.png')
        plt.savefig(savepath[0]+'_transposonplot_chrom'+chrom+'.png', dpi=400)
        plt.close()
    elif savefig == True and variable == "reads":
        savepath = os.path.splitext(bed_file)
        print('saving figure at %s' % savepath[0]+'_readplot_chrom'+chrom+'.png')
        plt.savefig(savepath[0]+'_readplot_chrom'+chrom+'.png', dpi=400)
        plt.close()
    else:
        plt.show()



#%%
if __name__ == '__main__':
    profile_plot(bed_file=bed_file, variable=variable, chrom=chrom, bar_width=bar_width, savefig=savefig)