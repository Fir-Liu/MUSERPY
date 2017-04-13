# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 11:32:17 2016

@author: Super
"""

import numpy as np
import misc as mi
#%%
import logging
reload(logging)
#import numpy as np
log_level = logging.DEBUG 
log_format = '%(levelname)s:%(message)s'
logger = logging.root
logger.basicConfig = logging.basicConfig(
#                                filename='myProgramLog.txt',
                                level=log_level, 
                                format=log_format)
#%%

def frstop(cordata,delay_sys,array='m1'):
    """
    Return the fringstopped xcor-data set. 
    
    xcordata: from rdin.rdraw()
    delay_par: which are extracted from raw-data frames, which is 
                the summation of system delay and gometric delay.
                shape=(Nant, Nts)
    delay_sys: measured system delay. shape=(Nant, Nts)
    array: 'm1' for muser1, 'm2' for muser2.
    """
    Nant = 40 if array=='m1' else 60
    bl2ord = mi.bl_tab(Nant)
    Nts = cordata['time'].shape[0]
    chan = np.arange(0,16)
    FLOrf=(1600+chan*25+12.5)/1000. # RF LO in GHz
    FLOif=(chan*25+12.5+50.)/1000.  # IF LO in GHz  
#    delay_par = 10+np.random.rand(40,Nts)
#    delay_sys = np.random.randint(5,10,size=(44,1))
#    delay_sys1 = np.repeat(delay_sys[0:40],Nts,axis=1)
    delay_geo = cordata['dly'] - delay_sys
    for i in range(Nant):
        for j in range(i+1,Nant):
            for v in range(16):
                for t in range(Nts):
    #                tg0=fix(delay[tt,j])-fix(delay[tt,i])
                    tg = delay_geo[j,t] - delay_geo[i,t]
                    tg0 = np.floor(delay_geo[j,t]) - np.floor(delay_geo[i,t])
                    phai = 2*np.pi*(FLOrf[v]*tg - FLOif[v]*tg0)
#                    logger.debug('phai is: %f',phai)
#                    logger.debug('dout[0,0,0] is: %f',cordata['x'][0,0,t].real)
                    cordata['x'][bl2ord[i,j],v,t] = cordata['x'][bl2ord[i,j],v,t]\
                                                    *np.exp(np.complex(0,-phai))
#                    logger.debug('dout[0,0,0] after fr is: %f',cordata['x'][0,0,t].real)                                                    
    return cordata                                                                                            
            
    