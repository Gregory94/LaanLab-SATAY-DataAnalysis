# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a test script file.
"""

import os, sys
import numpy
import pysam

path = os.path.join('/home', 'gregoryvanbeek', 'Documents', 'data_processing')
filename = os.path.join('E-MTAB-4885.WT2.bam')

file = os.path.join(path,filename)
print('Running path: ', file)

if os.path.isfile(file):
    print('File exists.')
elif os.path.exists(file):
    print('File does not exist, but path does exists.')
else:
    print('Path does not exist.')



bamfile = pysam.AlignmentFile(file, 'rb')


for read in bamfile.fetch(bamfile.get_reference_name(0), 1, 30):
    print(read)

str(read).split('\t')