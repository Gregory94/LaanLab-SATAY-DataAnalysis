   #CREATE A DICTIONARY WITH ALL GENES (BOTH THE COMMON NAME AND THEIR SYSTEMATIC NAME) AND SAVE THEM WITH THEIR RESPECTIVE LENGHTS (IN TERMS OF BP WHICH IS DEFINED AS bp=aa*3)

def list_gene_names(gene_information_file = None):
    '''Create a list of all known gene names and their aliases as listed on SGD (or as provided as an optional input file)
    Input is a standard file downloaded from https://www.uniprot.org/docs/yeast.
    Output is list of all genes, which also includes all the aliases (if they exists).
    '''

    if gene_information_file == None:
        gene_information_file = r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt'

    gene_name_list = []
    gene_counter = 0
    with open(gene_information_file) as f:
        lines = f.readlines()
        for i in range(58,len(lines)-6):    #THE GENES START AT LINE 58 AND STOP 6 LINES BEFORE THE END OF THE FILE.
            n=0
            l = lines[i]

            extra_columns = l.count(';')    #COUNT HOW MANY TIMES ';' OCCURS IN A LINE. THIS IS NEEDED TO GET THE RIGHT COLUMNS AS SOMETIMES ALIASES OF GENES ARE PRESENTED IN EXTRA COLUMNS
            l_short = ' '.join(l.split())
            l_list = l_short.split(' ')
            
            gene_name_list.append(l_list[0].strip(';'))
            gene_name_list.append(l_list[1+extra_columns].strip(';'))

            if l_list[1+extra_columns] == 'GAG' or l_list[1+extra_columns] == 'POL':    #THESE ARE SEQUENCES THAT SOMETIMES OCCUR WHICH HAVE TO BE IGNORED.
                extra_columns = extra_columns + 1
            if extra_columns > 0:
                for n in range(extra_columns):
                    gene_name_list.append(l_list[1+n].strip(';'))
            gene_counter += 1

    print('Number of genes found in file = ',gene_counter)
    return(gene_name_list)

if __name__ == '__main__':
    list_gene_names()