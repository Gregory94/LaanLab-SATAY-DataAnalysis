# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 09:38:16 2021

@author: gregoryvanbeek
"""

import sys

def samflags(flag=0, verbose=True):
    '''
    This script converts a decimal flag to binary and get the corresponding properties according to the sam-flag standard.
    The code is based on the explanation given here https://davetang.org/muse/2014/03/06/understanding-bam-flags/
    For manual checking sam flags, check https://broadinstitute.github.io/picard/explain-flags.html
    The input is a decimal number.
    The flags can be used from both sam and bam files (those should be the same).
    '''


    flag_binary = format(flag, '012b') # '#012b' to start the string with '0b'. 12 indicated that the string has length 12.

    if verbose==True:
        print('Entered decimal flag = %i' % flag)
        print('Corresponding binary flag = %s' % flag_binary)
        print('')


    if flag == 0:
        flagprop_list = ['read forward strand']

    else: #IF ANY LINE IN PROP_DICT IS CHANGED, ALSO CHANGE THIS IN THE TRANSPOSONMAPPING_SATAY.PY SCRIPT.
        prop_dict = {1: 'read paired',
                     2: 'read mapped in proper pair',
                     3: 'read unmapped',
                     4: 'mate unmapped',
                     5: 'read reverse strand',
                     6: 'mate reverse strand',
                     7: 'first in pair',
                     8: 'second in pair',
                     9: 'not primary alignment',
                     10: 'read fails platform/vendor quality checks',
                     11: 'read is PCR or optical duplicate',
                     12: 'supplementary alignment'}
    
    
        counter = 1
        flagprop_list = []
        for b in reversed(flag_binary):
            if int(b) == 1:
                flagprop_list.append(prop_dict.get(counter))
            counter += 1


    if verbose==True:
        print('PROPERTIES:')
        if not flagprop_list == []:
            for prop in flagprop_list:
                print(prop)
            del prop
        else:
            print('no properties found')


    return(flag_binary, flagprop_list)




if __name__ == '__main__':
    flag_binary, flagproperties = samflags(int(sys.argv[1]))
