# -*- coding: utf-8 -*-
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
def transposon_profile(chrom=None, bar_width=None, bed_file=None):
    '''This function creates a bar plot along the entire genome.
    The height of each bar represents the number of transposons at the genomic position indicated on the x-axis.
    The input is as follows: which chromosome (indicated by roman numeral), bar_width, bed_file.
    The bar_width determines how many basepairs are put in one bin. Little basepairs per bin may be slow. Too many basepairs in one bin and possible low transposon areas might be obscured.
    The bed_file is one of the files created by the Matlab code from the kornmann-lab.
    The background of the graph is color coded to indicate areas that code for genes.
    For this a list for essential genes is needed (used in 'list_known_essentials' function) and a .gff file is required (for the functions in 'chromosome_and_gene_positions.py') and a list for gene aliases (used in the function 'gene_aliases')
    '''


#%%
    bed_file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed"
    bar_width=None

    gff_file = os.path.join(file_dirname,'Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
    essential_genes_files = [os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_1.txt'),
                            os.path.join(file_dirname,'Data_Files','Cerevisiae_EssentialGenes_List_2.txt')]


#    counter = 0
#    for chrom in ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'XI', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']:

    chrom_list = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']
    
    chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)
    
    
    summed_chr_length_dict = {}
    summed_chr_length = 0
    for c in chrom_list:
        summed_chr_length_dict[c] = summed_chr_length
        summed_chr_length += chr_length_dict.get(c)    

    
    l_genome = 0
    for chrom in chrom_list:
        l_genome += int(chr_length_dict.get(chrom))
    print('Genome length: ',l_genome)
    if bar_width == None:
        bar_width = l_genome/1000
    
    

    gene_pos_dict = gene_position(gff_file)
    genes_currentchrom_pos_list = [k for k, v in gene_pos_dict.items()]
    genes_essential_list = list_known_essentials(essential_genes_files)


    with open(bed_file) as f:
        lines = f.readlines()
    

    chrom_names_dict, chrom_start_index_dict, chrom_end_index_dict= chromosome_name_bedfile(lines)

    alltransposoncounts_list = np.zeros(l_genome)
    for line in lines[chrom_start_index_dict.get("I"):chrom_end_index_dict.get("XVI")+1]:
        line = line.strip('\n').split()
        chrom_name = [k for k,v in chrom_names_dict.items() if v == line[0].replace("chr",'')][0]
        alltransposoncounts_list[summed_chr_length_dict.get(chrom_name) + int(line[1])] += 1


    alltransposoncounts_binnedlist = []
    val_counter = 0
    sum_values = 0
    for n in range(len(alltransposoncounts_list)):
        if int(val_counter % bar_width) != 0:
            sum_values += alltransposoncounts_list[n]
        elif int(val_counter % bar_width) == 0:
            alltransposoncounts_binnedlist.append(sum_values)
            sum_values = 0
        val_counter += 1
    alltransposoncounts_binnedlist.append(sum_values)
    
    allinsertionsites_list = np.linspace(0,l_genome,int(l_genome/bar_width+1))




    plt.figure(figsize=(17.0,6.0))
    grid = plt.GridSpec(20, 1, wspace=0.0, hspace=0.0)

    textsize = 20
    textcolor = "#000000"
    binsize = bar_width
    ax = plt.subplot(grid[0:19,0])
#    for gene in genes_currentchrom_pos_list:
#        if not gene_pos_dict.get(gene)[0] == 'Mito':
#            gene_start_pos = summed_chr_length_dict.get(gene_pos_dict.get(gene)[0]) + int(gene_pos_dict.get(gene)[1])
#            gene_end_pos = summed_chr_length_dict.get(gene_pos_dict.get(gene)[0]) + int(gene_pos_dict.get(gene)[2])
#            if gene in genes_essential_list:
#                ax.axvspan(gene_start_pos,gene_end_pos,facecolor="#BBE6AA",alpha=0.8)
#            else:
#                ax.axvspan(gene_start_pos,gene_end_pos,facecolor="#F6A089",alpha=0.8)
    ax.bar(allinsertionsites_list,alltransposoncounts_binnedlist,width=binsize,color="#00918f")
    ax.grid(False)
    ax.set_xlim(0,l_genome)

    for chrom in summed_chr_length_dict:
        ax.axvline(x = summed_chr_length_dict.get(chrom), linestyle='-', color=(0.9,0.9,0.9,1.0))

    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False) # labels along the bottom edge are off
    plt.ylabel('Transposon Count', fontsize=textsize, color=textcolor)
    plt.yticks([], [])
#    ax.axis("off")

    axc = plt.subplot(grid[19,0])
    for gene in genes_currentchrom_pos_list:
        if not gene_pos_dict.get(gene)[0] == 'Mito':
            gene_start_pos = summed_chr_length_dict.get(gene_pos_dict.get(gene)[0]) + int(gene_pos_dict.get(gene)[1])
            gene_end_pos = summed_chr_length_dict.get(gene_pos_dict.get(gene)[0]) + int(gene_pos_dict.get(gene)[2])
            if gene in genes_essential_list:
                axc.axvspan(gene_start_pos,gene_end_pos,facecolor="#00F28E",alpha=0.8)
            else:
                axc.axvspan(gene_start_pos,gene_end_pos,facecolor="#F20064",alpha=0.8)
    axc.set_xlim(0,l_genome)
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
if __name__ == '__main__':
    transposon_profile(chrom='XV', bar_width=2000, bed_file=r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder_S288C\align_out\ERR1533147_trimmed.sorted.bam.bed")