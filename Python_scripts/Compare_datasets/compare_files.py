# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:25:42 2020

@author: gregoryvanbeek
"""

import os
import warnings


wigfile1 = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533147_trimmed.sorted.bam.wig"
wigfile2 = r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\analysis_matlab_kornmanncode\ERR1533147_trimmed.sorted.bam.wig"

def compare_wigfiles(wigfile1, wigfile2):
    '''
    '''
    
    assert os.path.isfile(wigfile1), "Wig file not found at location %s" % wigfile1
    assert os.path.isfile(wigfile2), "Wig file not found at location %s" % wigfile2
    
    with open(wigfile1) as f:
        lines1 = f.readlines()

    with open(wigfile2) as f:
        lines2 = f.readlines()

    
    print(len(lines1))
    print(len(lines2))
    
    if len(lines1) < len(lines2):
        warnings.warn("Length of %s is smaller than length of %s" % (wigfile1, wigfile2))
    elif len(lines1) > len(lines2):
        warnings.warn("Length of %s is smaller than length of %s" % (wigfile2, wigfile1))


    line1_dict = {}
    line2_dict = {}
    for i in range(1,len(lines1)):
        if lines1[i].lower().startswith('variable'):
            chrom1 = lines1[1].replace('\n','')[-10:-1]
        else:
            line1 = lines1[i].replace('\n','').split(' ')
            line1_dict[chrom1+'_'+line1[0]] = int(line1[1])

    for i in range(1,len(lines2)):
        if lines2[i].lower().startswith('variable'):
            chrom2 = lines2[1].replace('\n','')[-10:-1]
        else:
            line2 = lines2[i].replace(' \n','').split(' ')
            line2_dict[chrom2+'_'+line2[0]] = int(line2[1])


    difference1_dict = {}
    for location in line1_dict:
        if not line1_dict.get(location) == line2_dict.get(location):
            difference1_dict[location] = [line1_dict.get(location), line2_dict.get(location)]

    difference2_dict = {}
    for location in line2_dict:
        if not line1_dict.get(location) == line2_dict.get(location):
            difference1_dict[location] = [line1_dict.get(location), line2_dict.get(location)]


#    print("Number of differences found: %i" % difference_counter)




if __name__ == '__main__':
    compare_wigfiles(r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533147_trimmed.sorted.bam.wig",
                     r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\analysis_matlab_kornmanncode\ERR1533147_trimmed.sorted.bam.wig")