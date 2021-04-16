'''
This script is based on gene_names.py
It creates a text file containing all protein names in Saccharomyces Cerevisiae to run in the geneSCF software.
The results are stored at the same location as the location of this script.

For running geneSCF, the following example command can be used when in the geneSCF folder:
./geneSCF -m=normal -i=/home/gregoryvanbeek/Documents/Software/geneSCF/S_Cerevisiae_protein_oln_name_full_genome.txt -o=/home/gregoryvanbeek/Documents/Software/geneSCF/ -db=GO_all -bg=6726 -p=no -org=sgd -t=sym
Use -m=update to update the gene database.

Use the following command to create a database:
./prepare_database -db=GO_all -org=sgd
'''

import os

def list_gene_names(gene_information_file = None):
    '''Create a list of all known gene names and their aliases as listed on SGD (or as provided as an optional input file)
    Input is a standard file downloaded from https://www.uniprot.org/docs/yeast.
    Output is list of all genes, which also includes all the aliases (if they exists).
    '''

    file_dirname = os.path.dirname(os.path.abspath('__file__'))
    if gene_information_file == None:
        gene_information_file = os.path.join(file_dirname,'..','..','data_files','Yeast_Protein_Names.txt')

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
    from datetime import date
    current_date = date.today()


    oln_saving_file = r"S_Cerevisiae_protein_oln_name_full_genome.txt"
    gene_oln_list_sorted = sorted(gene_oln_list)
    with open(oln_saving_file, 'w') as f:
#        f.write("org=S. Cerevisiae ; type=Genomic ; naming='oln' ; source file='Yeast_Protein_Names.txt' ; creation date=%s using 'gene_names.py'\n" % current_date)
        for oln_name in gene_oln_list_sorted:
            f.write("%s\n" % oln_name)



#    sgd_saving_file = r"S_Cerevisiae_protein_designation_name_full_genome.txt"
#    gene_sgd_list_sorted = sorted(gene_sgd_list)
#    with open(sgd_saving_file, 'w') as f:
#        
#        f.write("org=S. Cerevisiae ; type=Genomic ; naming='designation' ; source file='Yeast_Protein_Names.txt' ; creation date=%s using 'gene_names.py'\n" % current_date)
#        for sgd_name in gene_sgd_list_sorted:
#            f.write("%s\n" % sgd_name)

###
    print('Number of genes found in file = ',gene_counter)
    return(gene_name_list)



if __name__ == '__main__':
    list_gene_names()
