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

def ind2bl(flind,nant):
    """
       Given flind from flagind() and Nant,
       return flagged baseline cross correlation pair like 
       [(1,2),(2,3)...]
    """
    bl2ord = bl_tab(nant)
    rclist = []
    for sk in list(flind):
        indfl = np.where(bl2ord==sk)
        rclist.append(zip(indfl[0],indfl[1])[0])    
    return rclist
    
def flagind(Nant,flaglist):
    """
        Given one antenna flaglist,        
        return flag index array of 1-d flattened baseline-ordered index array
        0-1:[0], 0-2:[1],0-3:[2],1-2:[3], 1-3:[4],2-3:[5]        
        If the flagged baselines are: 0-1, 1-3, the flagind list is: [0 4] 
        
    """
#    Nant = 40
    bl2ord = bl_tab(Nant)
    s1 = np.zeros((len(flaglist),Nant))
    for idx,ind in enumerate(flaglist):
        s1[idx,:] = bl2ord[ind,:]
    s = s1[0,:]    
    for j in range(s1.shape[0]-1):
        s = np.array(list(set(s)|set(s1[j+1,:])))
    sind = np.array([k for k in s if k!=-1],dtype='u2')
    return sind



frdict1 = {'0.4G-0.8G':0,'0.8G-1.2G':1,'1.2G-1.6G':2, \
          '1.6G-2.0G':3}
frdict2 = {'2.0G-2.4G':0,'2.4G-2.8G':1, \
          '2.8G-3.2G':2,'3.2G-3.6G':3, '3.6G-4.0G':4,\
          '4.0G-4.4G':5, \
          '4.4G-4.8G':6,'4.8G-5.2G':7,'5.2G-5.6G':8, \
          '5.6G-6.0G':9,'6.0G-6.4G':10,'6.4G-6.8G':11, \
          '6.8G-7.2G':12,'7.2G-7.6G':13,'7.6G-8.0G':14, \
          '8.0G-8.4G':15,'8.4G-8.8G':16,'8.8G-9.2G':17, \
          '9.2G-9.6G':18,'9.6G-10.0G':19,'10.0G-10.4G':20, \
          '10.4G-10.8G':21,'10.8G-11.2G':22,'11.2G-11.6G':23, \
          '11.6G-12.0G':24,'12.0G-12.4G':25,'12.4G-12.8G':26, \
          '12.8G-13.2G':27,'13.2G-13.6G':28,'13.6G-14.0G':29,\
          '14.0G-14.4G':30,'14.4G-14.8G':31,'14.6G-15.0G':32}

infrdict1 = {v:k for k,v in frdict1.iteritems()}
infrdict2 = {v:k for k,v in frdict2.iteritems()}

def fr2ord(rffreq,array='m1'):
    ''' Returns the value in list range(0,33) that is responding to
        the frequency band like '0.4G-0.8G'...
    '''
    ord = frdict1[rffreq] if array=='m1' else frdict2[rffreq]
    return ord
              
def ord2fr(rfind,array='m1'):
    ''' Returns the frequency band like '0.4G-0.8G'...
        responding to list range(0,33) 
    '''
    rffreq = infrdict1[rfind] if array=='m1' else infrdict2[rfind]
    return rffreq    