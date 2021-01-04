# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 09:38:16 2021

@author: gregoryvanbeek
"""

def samflag(flag=0):
    '''
    '''

    flag_binary = format(flag, '012b') # '#012b' to start the string with '0b'. 12 indicated that the string has length 12.

    print('Entered decimal flag = %i' % flag)
    print('Corresponding binary flag = %s' % flag_binary)
    print('')


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
#    prop1 = 'read paired'
#    prop2 = 'read mapped in proper pair'
#    prop3 = 'read unmapped'
#    prop4 = 'mate unmapped'
#    prop5 = 'read reverse strand'
#    prop6 = 'mate reverse strand'
#    prop7 = 'first in pair'
#    prop8 = 'second in pair'
#    prop9 = 'not primary alignment'
#    prop10 = 'read fails platform/vendor quality checks'
#    prop11 = 'read is PCR or optical duplicate'
#    prop12 = 'supplementary alignment'

    counter = 1
    flagprop_list = []
    for b in reversed(flag_binary):
        if int(b) == 1:
            flagprop_list.append(prop_dict.get(counter))
        counter += 1
        
        return(flag_binary, flagprop_list)




if __name__ == '__main__':
    flag_binary, flagproperties = samflag(flag=2)
