# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:18:54 2016

@author: Super
"""
import numpy as np

def bl_tab(nant):
    ''' Returns a two-dimensional array bl2ord that will translate
        a pair of antenna indexes (antenna number - 1) to the ordinal
        number of the baseline in the 'x' key.  Note bl2ord(i,j) = bl2ord(j,i),
        and bl2ord(i,i) = -1.
    '''
    bl2ord = np.ones((nant,nant),dtype='int')*(-1)
    k = 0
    for i in range(nant-1):
        for j in range(i+1,nant):
            bl2ord[i,j] = k
            bl2ord[j,i] = k
            k+=1
    return bl2ord