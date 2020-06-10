   #CREATE A DICTIONARY WITH ALL GENES (BOTH THE COMMON NAME AND THEIR SYSTEMATIC NAME) AND SAVE THEM WITH THEIR RESPECTIVE LENGHTS (IN TERMS OF BP WHICH IS DEFINED AS bp=aa*3)

def gene_length_bp(gene_information_file = None):
    '''Create a dictionary for the gene length of all genes, in terms of basepairs, as listed on SGD.
    Input is a standard file downloaded from https://www.uniprot.org/docs/yeast.
    The length of the genes is determined by multiplying the number of aminoacids for the corresponding proteins by 3.
    Output is a dictionary with all genes (including aliases) as keys with values representing their basepair length.
    '''

    if gene_information_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        gene_information_file = os.path.join(file_dirname,'..','Data_Files','Yeast_Protein_Names.txt')

    gene_length_dict = {}
    with open(gene_information_file) as f:
        lines = f.readlines()
        for i in range(58,len(lines)-6):    #THE GENES START AT LINE 58 AND STOP 6 LINES BEFORE THE END OF THE FILE.
            n=0
            l = lines[i]

            extra_columns = l.count(';')    #COUNT HOW MANY TIMES ';' OCCURS IN A LINE. THIS IS NEEDED TO GET THE RIGHT COLUMNS AS SOMETIMES ALIASES OF GENES ARE PRESENTED IN EXTRA COLUMNS
            l_short = ' '.join(l.split())
            l_list = l_short.split(' ')
            
            if l_list[1+extra_columns] == 'GAG' or l_list[1+extra_columns] == 'POL':    #THESE ARE SEQUENCES THAT SOMETIMES OCCUR WHICH HAVE TO BE IGNORED.
                extra_columns = extra_columns + 1
            gene_length_dict[l_list[0].strip(';')] = int(l_list[5+extra_columns])*3 #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE NAME (e.g. Cdc42)
            gene_length_dict[l_list[1+extra_columns]] = int(l_list[5+extra_columns])*3  #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE SYSTEMATIC NAME (e.g. YLR229C)
            if extra_columns > 0:
                for n in range(extra_columns+1):
                    gene_length_dict[l_list[0+n].strip(';')] = int(l_list[5+extra_columns])*3  #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE ALIASES (IF PRESENT).

    return(gene_length_dict)


def gene_length_aa(gene_information_file = None):
    '''Create a dictionary for the gene length of all genes, in terms of amino acids, as listed on SGD.
    Input is a standard file downloaded from https://www.uniprot.org/docs/yeast.
    The length of the genes is determined by multiplying the number of aminoacids for the corresponding proteins by 3.
    Output is a dictionary with all genes (including aliases) as keys with values representing their amino acid length.
    '''

    if gene_information_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        gene_information_file = os.path.join(file_dirname,'..','Data_Files','Yeast_Protein_Names.txt')

    gene_length_dict = {}
    with open(gene_information_file) as f:
        lines = f.readlines()
        for i in range(58,len(lines)-6):    #THE GENES START AT LINE 58 AND STOP 6 LINES BEFORE THE END OF THE FILE.
            n=0
            l = lines[i]

            extra_columns = l.count(';')    #COUNT HOW MANY TIMES ';' OCCURS IN A LINE. THIS IS NEEDED TO GET THE RIGHT COLUMNS AS SOMETIMES ALIASES OF GENES ARE PRESENTED IN EXTRA COLUMNS
            l_short = ' '.join(l.split())
            l_list = l_short.split(' ')
            
            if l_list[1+extra_columns] == 'GAG' or l_list[1+extra_columns] == 'POL':    #THESE ARE SEQUENCES THAT SOMETIMES OCCUR WHICH HAVE TO BE IGNORED.
                extra_columns = extra_columns + 1
            gene_length_dict[l_list[0].strip(';')] = int(l_list[5+extra_columns]) #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE NAME (e.g. Cdc42)
            gene_length_dict[l_list[1+extra_columns]] = int(l_list[5+extra_columns])  #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE SYSTEMATIC NAME (e.g. YLR229C)
            if extra_columns > 0:
                for n in range(extra_columns+1):
                    gene_length_dict[l_list[0+n].strip(';')] = int(l_list[5+extra_columns])  #DETERMINE GENE LENGTH AND STORE WITH THE CORRESPONDING GENE ALIASES (IF PRESENT).

    return(gene_length_dict)

if __name__ == '__main__':
    gene_length_aa()