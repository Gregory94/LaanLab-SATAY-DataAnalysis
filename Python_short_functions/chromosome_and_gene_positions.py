
def chromosome_position(gff_file = None):
    '''Get the start and end position of each chromosome and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    Output are three dictionaries for length, start and end position. All dictionaries have keys representing the chromosome number in roman numerals.
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

    return('chr_length_dict':chr_length_dict, 'chr_start_pos_dict':chr_start_pos_dict, 'chr_end_pos_dict':chr_end_pos_dict)




def gene_position(gff_file = None):
    '''Get the start and end position of each gene and determine their respective length.
    Input is a .gff file downloaded from https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index
    '''

    if gff_file == None:
        gff_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Saccharomyces_cerevisiae.R64-1-1.99.gff3'


    with open(gff_file) as f:
        lines = f.readlines()





if __name__ == '__main__':
    chromosome_position()
