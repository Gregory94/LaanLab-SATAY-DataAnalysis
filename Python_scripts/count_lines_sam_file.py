# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:53:23 2020

@author: gregoryvanbeek
"""

import os

def linecounter_sam(samfile, text):
    
    assert os.path.isfile(samfile), "Samfile not found at location %s" % samfile

    line_counter = 0   
    with open(samfile) as f:
        for line in f:
            if text in line:
                line_counter += 1

    print(line_counter)

if __name__ == '__main__':
    linecounter_sam(r"C:\Users\gregoryvanbeek\Documents\testing_site\wt1_testfolder\align_out\ERR1533147_trimmed.sam", '0\tref|NC_001133|\t107')