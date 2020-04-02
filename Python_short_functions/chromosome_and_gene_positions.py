'''This module includes two functions, 'chromosome_position' and 'gene_position'.
Both do not require an input with the full path to the file 'Saccharomyces_cerevisiae.R64-1-1.99.gff3' (which can be downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index).
If no input is given, the program automatically searches for the file at r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3'.
'''
def chromosome_position(gff_file = None):
    '''Get the start and end position of each chromosome and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    Output are three dictionaries for length, start and end position. All dictionaries have keys representing the chromosome number in roman numerals.
    To get all dictionaries, use: 'a,b,c, chromosome_and_gene_position.chromosome_position()'.
    'a' = chromosome length
    'b' = chromosome start position
    'c' = chromosome end position
    '''
    if gff_file == None:
        gff_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3'

    #CONVERT ROMAN NUMERALS TO ARABIC NUMERALS
    roman_to_arabic_dict = {}
    roman_nums_list = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
    arabic_counter = 1
    for roman in roman_nums_list:
        roman_to_arabic_dict[roman] = arabic_counter
        arabic_counter += 1

    #GET END POSITIONS OF THE CHROMOSOMES FROM THE GFF FILE AND STORE THEM IN A DICTIONARY
    chr_length_dict = {}
    chr_length_list = []
    with open(gff_file) as f:
        line_counter = 0
        next(f)
        while line_counter < 17:
            lines = f.readline()
            chr_line_list = lines.strip('\n').replace(' ','\t').split('\t')
            chr_number = chr_line_list[3]
            if chr_number != 'Mito':
                chr_length = int(chr_line_list[5])
                chr_length_list.append(chr_length)
                chr_length_dict[chr_number] = chr_length
            line_counter += 1
    
    #DETERMINE START AND END POSITION OF EACH OF THE CHROMOSOMES
    chr_start_pos_dict = {}
    chr_end_pos_dict = {}
    counter = 0
    for roman in roman_nums_list:
        chr_start_pos_dict[roman] =  sum(chr_length_list[:counter])+1
        chr_end_pos_dict[roman] = sum(chr_length_list[:counter+1])
        if roman == 'I':
            chr_start_pos_dict[roman] = 1
        counter += 1

    return(chr_length_dict, chr_start_pos_dict, chr_end_pos_dict)




def gene_position(gff_file = None):
    '''Get the start and end position of each gene and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    Output is a dictionary that includes all gene names as keys. The values are lists with two numbers.
    The first is the chromosome number the gene belong to, the second is the start position, the third is the end position of the gene in terms of basepairs.
    '''

    if gff_file == None:
        gff_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3'

    gene_pos_dict = {}
    with open(gff_file) as f:
        for line in f:
            line_list = line.split('\t')
            if len(line_list) > 2:
                if line_list[2] == 'gene':
                    gene_chr = line_list[0]
                    gene_start = line_list[3]
                    gene_end = line_list[4]
                    gene_position = [gene_chr,gene_start,gene_end]
                    
                    gene_name_string = line_list[8].split(';')[0]
                    gene_name = gene_name_string.split(':')[1]
                    
                    gene_pos_dict[gene_name] = gene_position

    return(gene_pos_dict)

if __name__ == '__main__':
    chromosome_position()
