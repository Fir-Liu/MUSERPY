# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 13:05:05 2016

@author: Super
"""

#%%        
import numpy as np
import misc as mc
import astropy.constants as con
from astropy import units as u
import matplotlib.pyplot as plt
import numpy.ma as ma
import img
import rdin as ri

#
arraycfgfile = r'..\config\LFA_Ant_pos.txt'
obs_freq = 1.7025e+9 
#site_lat = 42.2 * np.pi / 180
#site_lon = 115.25 * np.pi / 180
site_hr = (4*15) * np.pi / 180
site_dec = 21.71*np.pi/180
antflaglist = [4,17,29]
# Set imaging parameters
Imdim = 512
uvdarr = np.zeros((Imdim,Imdim))
Imsize = 1*u.deg      
uvres = np.floor((1/Imsize.to(u.rad))).value
rowt = 4*uvres


def imsizecal(FOV,UVMAX,NBP=3):
    """
    FOV: Field of View, or Image physical size in degree.
    UVMAX: absolute uvspan maximum
    NBP: pixel points to representing synthesized beam.
          NBP >=2
    
    Return image size M, so image dim is MxM
    """
    return np.pi*FOV/180*UVMAX*NBP

#%%    
# need to flag antenna and make uvwbl conj symmetry
uvwbl,flbllist = img.muvw('m1',antflaglist,obs_freq,site_hr,site_dec)

FOV = 1
UVMAX=max(uvwbl['u'].max(),uvwbl['v'].max())

M = int(np.ceil(imsizecal(FOV,UVMAX)))
#%%
# read correlation data in
filename = r'C:\MUSER_DATA\CSRH_20141111-122131_98335183'
trange = ['2014-11-11 12:22:35', '2014-11-11 12:22:45']
dout = ri.rdraw(filename,trange,tdec=0,array='m1', 
          nant=40,rfind=0,pol='LL')
#%%
# Flag uvdata
dout_flag = img.mflag('m1',uvwbl,dout['x'][:,0,0],antflaglist)          
#%%
# plot flagged uv distribution
u_flag = ma.getdata(uvwbl['u'][uvwbl['u'].mask])
v_flag = ma.getdata(uvwbl['v'][uvwbl['v'].mask])
plt.figure(3)
plt.scatter(uvwbl['u'],uvwbl['v'],c='b',marker='.')
plt.scatter(u_flag,v_flag,c='r',marker='s')
#plt.scatter()
plt.grid(True)
#%%
gdpts = np.ceil(max(uvwbl['u'].max(),uvwbl['v'].max())/uvres)+1
uvgspan = np.arange(-gdpts*uvres,(gdpts+1)*uvres,uvres)
M = 2*gdpts+1

#%%
# Set xcor dataset with numpy structured masked array
dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]
Np = len(ma.compressed(uvwbl['u']))
psf = np.array(np.zeros(Np),  dtype=dt)
psf['xk'] = np.ones(Np)
psf['wgtk'] = np.ones(Np)
psf['uk'] = ma.compressed(uvwbl['u'])
psf['vk']= ma.compressed(uvwbl['v'])
psf['wk']= ma.compressed(uvwbl['w'])

#%%
# Improving weighting performance

psf_wt = img.mweight(psf,rowt)    
#%%
# Gridding with discrete convolution
kervspan = 3*uvres
rms = uvres/2
kerv = ('g',rms) # gaussian function is kernel
xcorg = img.mgrid(psf_wt,uvgspan,uvres,kerv,kervspan)    
#%%    
Nf = int(M)
Nz = Imdim**2 - Nf**2
Nzl,Nzr = Nz/2+1, Nz/2
xcorg_ex = np.hstack((np.zeros(Nzl),xcorg['xg'],np.zeros(Nzr)))
xcorg_mat = xcorg_ex.reshape(Imdim,Imdim)
beam = np.fft.ifft2(np.fft.fftshift(xcorg_mat))

beam = np.fft.fftshift(np.real(beam))
plt.figure(2)
imbeam = plt.imshow(beam)
plt.colorbar()


#%%
%timeit -n1 -r1 foo(psf,rowt)    
%timeit -n1 -r1 gridconv(psf_wt,uvgspan,kerv,kervspan)    

#%%
## Weighting---critical part
# computes density of given region
#%timeit -n1 -r1 foo(mDUVW_wgt,Nbl)
#def foo(psf,rowt):
#    for k in range(psf.shape[0]/2):
#        for j in range(0,k)+range(k+1,psf.shape[0]):
#            if (np.abs(psf['uk'][j]-psf['uk'][k])<=rowt
#                and np.abs(psf['vk'][j]-psf['vk'][k])<=rowt):
#               psf['wgtk'][k] = psf['wgtk'][k] + 1 
#               psf['wgtk'][k+Nbl] = psf['wgtk'][k+Nbl] + 1            
#            else:       
#                pass       
#    return psf
##    return mDUVW_wgt
#antxyz = np.zeros((40,3))
#f=open(arraycfgfile,'r')
#header1 = f.readline()
#Nant=40
#Nbl = Nant*(Nant-1)/2
##k=0
#for ik,line in enumerate(f):
#    line = line.strip()
#    cline = np.array(line.split())
#    antxyz[ik,:] = cline.astype(np.float)
##    k=k+1
#f.close()
#%%
# Set location and observing parameters
# and compute (u,v,w)

#DXYZ = np.zeros((3,Nbl))
#bl2ord = mc.bl_tab(Nant)
#
#for j in range(Nant-1):
#    for i in range(j+1,Nant):
#        DXYZ[:,bl2ord[j,i]] = antxyz[j,:]-antxyz[i,:]
#nu2uv_ma = mc.nu2uv(site_lat, site_lon, site_hr, site_dec)
#DUVW = np.dot(nu2uv_ma, DXYZ)*obs_freq/con.c.value
#
#dtuvw = [('u',float),('v',float),('w',float)]
#uvwbl = np.zeros(Nbl,dtype=dtuvw)
#uvwbl['u'] = [ubl for ubl in DUVW[0,:]]
#uvwbl['v'] = [vbl for vbl in DUVW[1,:]]
#uvwbl['w'] = [wbl for wbl in DUVW[2,:]]
