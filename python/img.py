# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 10:33:02 2016

@author: Super
"""
import numpy as np
import misc as mi
import coord as cor
import astropy.constants as con
import numpy.ma as ma
import rdin as ri
#%%
# printing for debug
import logging
import imp
imp.reload(logging)
#import numpy as np
log_level = logging.DEBUG 
log_format = '%(levelname)s:%(message)s'
logger = logging.root
logger.basicConfig = logging.basicConfig(
#                                filename='myProgramLog.txt',
                                level=log_level, 
                                format=log_format)
logging.disable(logging.DEBUG)
#%%
def muvw(arrayname,flaglist,obs_freq,site_hr,site_dec,site_lat=42.2*np.pi/180, 
         site_lon=115.25 * np.pi / 180,configpath=r'..\config' ):
    """
        Read array configuration and compute (u,v,w)
    """
    m1file = 'LFA_Ant_pos.txt'
    m2file = 'HFA_Ant_pos.txt'
    configfile = configpath+'\\'+ m1file if arrayname == 'm1' else configpath+m2file
    Nant = 40 if arrayname == 'm1' else 60
    Nbl = Nant*(Nant-1)/2
    antxyz = np.zeros((Nant,3))
    f=open(configfile,'r')
    header1 = f.readline() # ignore the header
    for ik,line in enumerate(f):
        line = line.strip()
        cline = np.array(line.split())
        antxyz[ik,:] = cline.astype(np.float)
    f.close() 

    DXYZ = np.zeros((3,Nbl))
    bl2ord = mi.bl_tab(Nant)
    
    for j in range(Nant-1):
        for i in range(j+1,Nant):
            DXYZ[:,bl2ord[j,i]] = antxyz[j,:]-antxyz[i,:]
    nu2uv_ma = cor.nu2uv(site_lat, site_lon, site_hr, site_dec)
    DUVW = np.dot(nu2uv_ma, DXYZ)*obs_freq/con.c.value
    
#    dtuvw = [('u',float),('v',float),('w',float)]
    uvwbl = ma.masked_array(np.zeros(2*Nbl,dtype=ri.dtuvw))
    uvwbl['u'] = np.hstack((DUVW[0,:],-DUVW[0,:]))
    uvwbl['v'] = np.hstack((DUVW[1,:],-DUVW[1,:]))
    uvwbl['w'] = np.hstack((DUVW[2,:],-DUVW[2,:]))
    
    flindls = mi.flagind(Nant,flaglist)
    flindls_ex = np.hstack((flindls,flindls+Nbl)) # flag (-u,-v)
    for ifl in flindls_ex:
        uvwbl[ifl] = ma.masked               
#    uvwbl_flag = np.array(zip(ma.compressed(uvwbl['u']),
#                        ma.compressed(uvwbl['v']),
#                        ma.compressed(uvwbl['w'])),
#                        dtype=ri.dtuvw)
    flrclist = mi.ind2bl(flindls,Nant)
    return uvwbl,flrclist

def mflag(arrayname,uvwbl,uvdata,flaglist):
    """
        Return flagged uvdata                
    """
    Nant = 40 if arrayname == 'm1' else 60
    Nbl = Nant*(Nant-1)/2
#    logger.debug('Nbl is %d',Nbl)
    dout_mask = ma.masked_array(zip(uvwbl['u'],uvwbl['v'],uvwbl['w'],
                            np.hstack((uvdata[0:Nbl],np.conj(uvdata[0:Nbl]))),
                            np.ones(Nbl*2)),dtype=ri.dtk)
    flagarr = mi.flagind(Nant,flaglist)
    flagarr = np.hstack((flagarr,flagarr+Nbl)) # flag (-u,-v)
    for ifl in flagarr:
        dout_mask[ifl] = ma.masked               
#    Np = ma.count(dout_mask.dtype.names[0])   
#    logger.debug('Np is %d',Np)    
#    logger.debug('dout shape is %d',dout.shape[0])
    dout_flag = np.array(zip(ma.compressed(dout_mask['u']),
                        ma.compressed(dout_mask['v']),
                        ma.compressed(dout_mask['w']),
                        ma.compressed(dout_mask['x']),
                        ma.compressed(dout_mask['wgt'])),dtype=ri.dtk)
    return dout_flag
     


def mpsf(uvwbl,uvgspan,uvres,
         kerv,kervspan,wtregion,weight,taper=None):
     """
        Return weighted point spread function                
     """
     return 0
   
def mweight(uvdata,wtregion,weight='U',taper=None): 
    """
        Return weighting Point-Spread-Function(PSF)
        uvdata is numpy array-like data with dtype = [('u',float),('v',float),('w',float),('x',float/complex),('wgt',float)]
        weight: 'u'-uniform, 'n'-natural
    """
#    uvdata = ma.compressed(uvdata) if ma.isMA(uvdata) else uvdata
#    fdname = uvdata.dtype.names
    if weight == 'U': 
        for idx,uvargs in enumerate(uvdata): # searching over half of psf is enough
            C1 = np.logical_and(np.abs(uvdata['u']-uvargs[0])<wtregion,uvdata['u']!=uvargs[0])
            C2 = np.logical_and(np.abs(uvdata['v']-uvargs[1])<wtregion,uvdata['v']!=uvargs[1])
            kind = np.where(np.logical_and(C1,C2))[0]
            uvdata['wgt'][idx] = len(kind) if len(kind)>0 else uvdata['wgt'][idx]
        uvdata['x'] = uvdata['x']/uvdata['wgt']
    else:
        pass    
    return uvdata    

def mgrid2mat(data,dim):
#    dt = data.dtype.names
    data_st = np.sort(data,order=['ug','vg'])
    datamat = np.zeros((dim,dim),dtype='complex')
    rcind_uv = [(i,j) for j in range(dim) for i in range(dim)[::-1]]
    for idx,(i,j) in enumerate(rcind_uv):
        datamat[i,j] = data_st['xg'][idx]                
    return datamat
    
def mgrid(uvdata,uvgspan,uvres,kerv,kervspan):
    def kfun(xseq, yseq,ktype):
        """
            return sequence: [f(xseq[0],yseq[0]),...]
        """
        if ktype[0]=='g':
            rms = ktype[1]
            kf = lambda x,y : np.exp(-x**2/(2*rms**2)-y**2/(2*rms**2))    
            return map(kf,xseq,yseq)    
        else:
            return 0
#    uvdata = ma.compressed(uvdata) if ma.isMA(uvdata) else uvdata    
#    fdname = uvdata.dtype.names
#    dtx = [('ug',float),('vg',float),('xg',complex)]
    xcorg = np.array([(ux,vx,0) for vx in uvgspan[::-1] for ux in uvgspan],dtype = ri.dtg)
    N = len(xcorg)
    for idg,(ug,vg,_) in enumerate(xcorg[0:(N-1)/2+1]):
        inda = np.where(np.logical_and(np.abs(uvdata['v']-vg)<kervspan,
                           np.abs(uvdata['u']-ug)<kervspan))[0]
        xcorg['xg'][idg] = uvdata['x'][inda].dot(kfun(uvdata['u'][inda]-ug, 
                            uvdata['v'][inda]-vg,kerv)) if len(inda)>0 else xcorg['xg'][idg]
        xcorg['xg'][N-1-idg] = xcorg['xg'][idg]
    return xcorg
 
def UVzeropad(s,Ndim):
    """
        Zeropadding a 2-d space frequency matrix s, which causes
        dims of s extended to Ndim. s must be arranged in a Hemitian way.
        s[i,j]=np.conj(s[N-1-i,N-1-j])
        Zeropadding of s interpolates zeros at freqs near nyquiste frequency.
    """
    Ns = s.shape[0]
    Nzlu = int(np.ceil((Ndim-Ns)/2.))
    Nzrd = int(np.floor((Ndim-Ns)/2.))
    sz = np.pad(s,[Nzlu,Nzrd],'constant',constant_values=0)
    return sz    