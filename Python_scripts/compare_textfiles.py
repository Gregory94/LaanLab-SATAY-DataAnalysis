# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 10:53:48 2020

@author: gregoryvanbeek

THIS CODE READS TWO COMMA DELIMITED TEXTS FILES AND COMPARES THE DATA WITH EACH OTHER AND INDICATES DIFFERENCES.
THE ORDER IN WHICH THE PROGRAM SHOWS THE DIFFERENCES IN THE SAME ORDER AS WHICH THE PATHS ARE GIVEN.
"""

import os

#%% compare variable: tncoordinates
paths = [r"C:\Users\gregoryvanbeek\Desktop\Python_matlab_differences\tncoordinatescopy_matlab.txt",
         r"C:\Users\gregoryvanbeek\Desktop\Python_matlab_differences\tncoordinatescopy_python.txt"]

for path in paths:
    if not os.path.exists(path):
        print('The following path does not exists:')
        print(path)


file0 = open(paths[0])
lines_file0 = file0.readlines()
file0.close()
print('Length ', os.path.basename(paths[0]), 'is %i' %int(len(lines_file0)))

file1 = open(paths[1])
lines_file1 = file1.readlines()
file1.close()
print('Length ', os.path.basename(paths[1]), 'is %i' %int(len(lines_file1)))
N = len(lines_file1)

differences_list = []
for ii in range(0,N):
    line0 = [int(x) for x in lines_file0[ii].split(',')]
    line1 = [int(x) for x in lines_file1[ii].split(',')]
    
    if not line0 == line1:
        differences_list.append([ii,line0,line1])
        if abs(line0[1] - line1[1]) >= 2:
            print('Differences in insertion location is bigger than 2 at location ', ii)


#%% compare variable: aa
paths = [r"C:\Users\gregoryvanbeek\Desktop\Python_matlab_differences\aa_matlab.txt",
         r"C:\Users\gregoryvanbeek\Desktop\Python_matlab_differences\aa_python.txt"]

for path in paths:
    if not os.path.exists(path):
        print('The following path does not exists:')
        print(path)


file0 = open(paths[0])
temp_file0 = file0.readlines()
file0.close()
lines_file0 = []
for line in temp_file0:
    lines_file0.append(int(line))
print('Length ', os.path.basename(paths[0]), 'is %i' %int(len(lines_file0)))

file1 = open(paths[1])
temp_file1 = file1.readlines()
file1.close()
lines_file1 = temp_file1[0].split(',')
for i in range(0,len(lines_file1)):
    lines_file1[i] = int(lines_file1[i]) + 1# +1 to account that python starts counting at 0 and matlab at 0.
print('Length ', os.path.basename(paths[1]), 'is %i' %int(len(lines_file1)))

N = len(lines_file1)

differences_list = []
for ii in range(0,N):
    if not lines_file0[ii] == lines_file1[ii]:
        differences_list.append([ii,lines_file0[ii],lines_file1[ii]])



