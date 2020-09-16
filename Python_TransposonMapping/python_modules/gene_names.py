'''
'''

def list_gene_names(gene_information_file = None):
    '''Create a list of all known gene names and their aliases as listed on SGD (or as provided as an optional input file)
    Input is a standard file downloaded from https://www.uniprot.org/docs/yeast.
    Output is list of all genes, which also includes all the aliases (if they exists).
    '''

    if gene_information_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        gene_information_file = os.path.join(file_dirname,'..','..','..','data_files','Yeast_Protein_Names.txt')

    gene_name_list = [] #INLCUDES ALL GENE NAMES AND POTENTIAL ALIASES
    gene_oln_list = [] #INCLUDE ONLY THE OLN NAMING CONVENTION
    gene_sgd_list = [] #INCLUDE THE FIRST DESIGNATION NAME THAT IS IN THE INPUT LIST
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
            gene_sgd_list.append(l_list[0].strip(';'))

            gene_oln = l_list[1+extra_columns].strip(';') #GET THE OLN NAME
            if gene_oln == 'GAG' or gene_oln == 'POL': #CHECK WHETHER THE OLN IS 'GAG' OR 'POL'. IF YES, TAKE THE NEXT COLUMN
                gene_name_list.append(l_list[2+extra_columns].strip(';'))
                gene_oln_list.append(l_list[2+extra_columns].strip(';'))
            else:
                gene_name_list.append(gene_oln)
                gene_oln_list.append(gene_oln)

            if l_list[1+extra_columns] == 'GAG' or l_list[1+extra_columns] == 'POL':    #THESE ARE SEQUENCES THAT SOMETIMES OCCUR WHICH HAVE TO BE IGNORED.
                extra_columns = extra_columns + 1
            if extra_columns > 0:
                for n in range(extra_columns):
                    gene_name = l_list[1+n].strip(';')
                    if not gene_name == 'GAG' and not gene_name == 'POL':
                        gene_name_list.append(gene_name)
            gene_counter += 1

###SAVING OLN LIST
#    from datetime import date
#    current_date = date.today()


#    oln_saving_file = r"S_Cerevisiae_protein_oln_name_full_genome.txt"
#    gene_oln_list_sorted = sorted(gene_oln_list)
#    with open(oln_saving_file, 'w') as f:
#        
##        f.write("org=S. Cerevisiae ; type=Genomic ; naming='oln' ; source file='Yeast_Protein_Names.txt' ; creation date=%s using 'gene_names.py'\n" % current_date)
#        for oln_name in gene_oln_list_sorted:
#            f.write("%s\n" % oln_name)



#    sgd_saving_file = r"S_Cerevisiae_protein_designation_name_full_genome.txt"
#    gene_sgd_list_sorted = sorted(gene_sgd_list)
#    with open(sgd_saving_file, 'w') as f:
#        
#        f.write("org=S. Cerevisiae ; type=Genomic ; naming='designation' ; source file='Yeast_Protein_Names.txt' ; creation date=%s using 'gene_names.py'\n" % current_date)
#        for sgd_name in gene_sgd_list_sorted:
#            f.write("%s \n" % sgd_name)

###
    print('Number of genes found in file = ',gene_counter)
    return(gene_name_list)









def gene_aliases(gene_information_file = None):
    '''Create three dictionaries containing aliases for genes
    Input is the path to 'Protein_Names.txt' file downloaded from https://www.uniprot.org/docs/yeast.
    If no input is given the file is automatically searched for at r'X:\tnw\BN\LL\Shared\Gregory\Gene_Database\Yeast_Protein_Names.txt'.
    Output is three dictionaries:
    aliases_designation_dict = gene aliases for common names (e.g. Bem1 and Sro1)
    aliases_sgd_dict = gene aliases for the search names in SGD (e.g. Bem1 and S000000404)
    aliases_swissprot_dict = gene aliases for Swiss Prot.
    The keys of the dictionaries are the systematic names of the genes (e.g. YBR200W for Bem1)
    
    Search through lists to get corresponding key:
        [key for key, val in aliases.items() if 'TFC3' in val]
    '''

    if gene_information_file == None:
        import os
        file_dirname = os.path.dirname(os.path.abspath('__file__'))
        gene_information_file = os.path.join(file_dirname,'..','Data_Files','Yeast_Protein_Names.txt')


    aliases_designation_dict = {}
    aliases_sgd_dict = {}
    aliases_swissprot_dict = {}
    
    
    with open(gene_information_file) as f:
        lines = f.readlines()
        
        for i in range(58,len(lines)-6):    #THE GENES START AT LINE 58 AND STOP 6 LINES BEFORE THE END OF THE FILE.
            l = lines[i]
            l_list = l.split()
            
             
            alias_counter = 0
            gene_designation_list = []
            for names in l_list:
                if names.endswith(';'):
                    gene_designation_list.append(names.strip(';'))
                    alias_counter += 1
            gene_designation_list.append(l_list[alias_counter])
            
            
            oln_name = l_list[alias_counter+1]
            sgd_name = l_list[alias_counter+4]
            swissprot_name = l_list[alias_counter+2]
            if oln_name == 'GAG' or oln_name == 'POL':
                oln_name = l_list[alias_counter+2]
                sgd_name = l_list[alias_counter+5]
                swissprot_name = l_list[alias_counter+3]

            aliases_designation_dict[oln_name] = gene_designation_list
            aliases_sgd_dict[oln_name] = sgd_name
            aliases_swissprot_dict[oln_name] = swissprot_name



    return(aliases_designation_dict, aliases_sgd_dict, aliases_swissprot_dict)

if __name__ == '__main__':
    list_gene_names()
