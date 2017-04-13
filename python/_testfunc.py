# -*- coding: utf-8 -*-
"""
Created on Tue Sep 06 11:00:43 2016

@author: Super
"""
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
#logging.disable(logging.DEBUG)

import numpy as np
import misc as mc
import astropy.constants as con
from astropy import units as u
import matplotlib.pyplot as plt
import numpy.ma as ma
import img
import rdin as ri
from scipy import signal
import deconv as dv
from numpy.fft import fft,ifft,fft2,ifft2,fftshift,ifftshift


tmap = np.arange(0,100).reshape((10,10))
tker = np.around(kernel,2)    
etker,etmap = dv.effregion(tker,tmap,5,5)

#test etker and etmap over all (i,j) ,tmap.shape
btst = np.ones(tmap.shape,dtype='bool')

for i in range(tmap.shape[0]):
    for j in range(tmap.shape[1]):
        etker,etmap = dv.effregion(tker,tmap,i,j)
        ixs=zip(*np.where(etmap==i*10+j))[0]
        btst[i,j]=(etker[ixs]==etker.max())&(etker.shape==etmap.shape)

btst.all()        
