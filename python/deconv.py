# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 14:01:15 2016

@author: Super

Deconvolution method
CLEAN
MEM
"""

import numpy as np
from scipy import signal
from astropy.modeling import models, fitting

def effregion(kernel,smap,*idxc):
    """
    Move kernel's center to (i,j).Kernel's span could be out of the 
    border of smap, return the overlapping region of kernel and smap.
    """
    ksx,ksy = kernel.shape
    Ms,Ns = smap.shape
    xa,xb = idxc[0]-np.floor_divide(ksx,2),idxc[0]+np.floor_divide(ksx,2)
    ya,yb = idxc[1]-np.floor_divide(ksy,2),idxc[1]+np.floor_divide(ksy,2)
    osxa = 0 if (xa<0) else xa
    osxb = Ms-1 if (xb>Ms-1) else xb
    osya = 0 if (ya<0) else ya
    osyb = Ns-1 if (yb>Ns-1) else yb

    okxa = ksx-(ksx+xa) if (xa<0) else 0
    okxb = ksx-1-(xb-Ms+1) if (xb>Ms-1) else ksx-1
    okya = ksy-(ksy+ya) if (ya<0) else 0
    okyb = ksy-1-(yb-Ns+1) if (yb>Ns-1) else ksy-1
#    smap1 = smap[osxa:osxb+1,osya:osyb+1]
#    kernel1 = kernel[okxa:okxb+1,okya:okyb+1]
#    return kernel1,smap1
    return (okxa,okxb,okya,okyb),(osxa,osxb,osya,osyb)



def hclean(dirtyim, psf, thresh, niter, lpgain, window=None):
    """
    Hogbom CLEAN alogrithm
    :dirtyim:
    :psf:
    :thresh:
    :niter:
    :lpgain:
    :window:        
        
    """
    clmodel = np.zeros_like(dirtyim)
    res = np.copy(dirtyim)
    for n in range(niter):
        [mx,my] = np.unravel_index(np.fabs(res).argmax(),res.shape)
        clmodel[mx,my] += lpgain*res[mx,my]
        okxy,osxy= effregion(psf,res,mx,my)
        res[osxy[0]:osxy[1]+1,osxy[2]:osxy[3]+1] -= psf[okxy[0]:okxy[1]+1,okxy[2]:okxy[3]+1]*lpgain*res[mx,my]
        if (np.fabs(res).max()< thresh):
            break

    xg,yg = np.arange(-3,4),np.arange(-3,4)
    xpg,ypg = np.meshgrid(xg,yg)
    
    g_init = models.Gaussian2D(amplitude=1., x_mean=0, x_stddev=1.,
                               y_mean=0, y_stddev=1.,theta=0)
    fitg = fitting.LevMarLSQFitter()
    g2d = fitg(g_init, xpg, ypg,psf)     
    # Using 2-d Gaussian beam as clean beam to fit psf
    clbeam = g2d(xpg,ypg)                      
    clmap = signal.fftconvolve(clmodel,clbeam,'same')
        
    return clmap,res

def clclean(dirtyim, psf, window, thresh, niter, lpgain):
    return 0
    
