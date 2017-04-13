# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 12:00:16 2016

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
#%%
import numpy as np
import misc as mc
import astropy.constants as con
from astropy import units as u
import matplotlib.pyplot as plt
import numpy.ma as ma
import img
import rdin as ri
from scipy import signal
from astropy.modeling import models, fitting
import deconv as dv
from numpy.fft import fft,ifft,fft2,ifft2,fftshift,ifftshift

#%%
def GetDirtyBeam():
    """
        Return the dirty beam of specified U-V distribution    
    """
#    arraycfgfile = r'..\config\LFA_Ant_pos.txt'
    obs_freq = 1.7025e+9 
    #site_lat = 42.2 * np.pi / 180
    #site_lon = 115.25 * np.pi / 180
    site_hr = (4*15) * np.pi / 180
    site_dec = 21.71*np.pi/180
    antflaglist = [11,12,13,24,25,26,37,38,39]
    # Set imaging parameters
#    Imdim = 512
    #    uvdarr = np.zeros((Imdim,Imdim))
    Imsize = 1*u.deg      
    uvres = np.floor((1/Imsize.to(u.rad))).value
    rowt = 4*uvres
    
    # need to flag antenna and make uvwbl conj symmetry
    uvwbl,flbllist = img.muvw('m1',antflaglist,obs_freq,site_hr,site_dec)
    #    u_flag = ma.getdata(uvwbl['u'][uvwbl['u'].mask])
    #    v_flag = ma.getdata(uvwbl['v'][uvwbl['v'].mask])
    
    gdpts = np.ceil(max(uvwbl['u'].max(),uvwbl['v'].max())/uvres)+1
    uvgspan = np.arange(-gdpts*uvres,(gdpts+1)*uvres,uvres)
    M = int(2*gdpts+1)
#    dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]
    Np = len(ma.compressed(uvwbl['u']))
    psf = np.array(np.zeros(Np),  dtype=ri.dtk)
    psf['x'] = np.ones(Np)
    psf['wgt'] = np.ones(Np)
    psf['u'] = ma.compressed(uvwbl['u'])
    psf['v']= ma.compressed(uvwbl['v'])
    psf['w']= ma.compressed(uvwbl['w'])
    
    psf_wt = img.mweight(psf,rowt)    
    kervspan = 3*uvres
    rms = uvres/2
    kerv = ('g',rms) # gaussian function is kernel
    psfg = img.mgrid(psf_wt,uvgspan,uvres,kerv,kervspan)    
    #   
#    dt1 = [('ug',float),('vg',float),('xg',float)]
    psf_mat = img.mgrid2mat(psfg,M)
#    beam = np.abs(ifft2((psf_mat)))
    return psf_mat


#%%
#def GetDirtyBeamA():
#    """
#        Return the dirty beam of specified U-V distribution    
#        beamwidth is 3 pixel points. 
#    """
#    arraycfgfile = r'..\config\LFA_Ant_pos.txt'
obs_freq = 1.7025e+9 
#site_lat = 42.2 * np.pi / 180
#site_lon = 115.25 * np.pi / 180
site_hr = (4*15) * np.pi / 180
site_dec = 21.71*np.pi/180
antflaglist = [11,12,13,24,25,26,37,38,39]
# Set imaging parameters
FOV = 1*u.deg   # Field of View or Imaging Size   
Npb = 3
NFT = 512
uvres = np.floor((1/FOV.to(u.rad))).value # resolution of ?

# need to flag antenna and make uvwbl conj symmetry
uvwbl,flbllist = img.muvw('m1',antflaglist,obs_freq,site_hr,site_dec)

gdpts = np.ceil(max(uvwbl['u'].max(),uvwbl['v'].max())/uvres)+1
umax = max(uvwbl['u'].max(),uvwbl['v'].max())

Npb1 = NFT/((FOV.to(u.rad)).value*umax)
ud = 2*umax/NFT
ugd = np.ceil(ud)
uvgspan = np.arange(-NFT/2*ugd,(NFT/2+1)*ugd,ugd)
rowt = 4*ugd
M = int(2*gdpts+1)
#    dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]
Np = len(ma.compressed(uvwbl['u']))
psf = np.array(np.zeros(Np),  dtype=ri.dtk)
psf['x'] = np.ones(Np)
psf['wgt'] = np.ones(Np)
psf['u'] = ma.compressed(uvwbl['u'])
psf['v']= ma.compressed(uvwbl['v'])
psf['w']= ma.compressed(uvwbl['w'])

psf_wt = img.mweight(psf,rowt)    
kervspan = 3*ugd
rms = ugd/2
kerv = ('g',rms) # gaussian function is kernel
psfg = img.mgrid(psf_wt,uvgspan,ugd,kerv,kervspan)    
#%%   
psf_mat = img.mgrid2mat(psfg,NFT+1)
psf_mat = psf_mat[0:NFT,0:NFT]
#    return psf_mat



    #%%    
def PointSource(imsize,region,*loc):
    """
        Return a N-point-source image      
    """
    image0 = np.zeros((imsize,imsize))
#    lent = len(loc)
#    logger.debug(list(loc))
    for v in list(loc):
#        logger.debug(v)
        image0[v[0]:v[0]+region+1,v[1]:v[1]+region+1] = 10 
    
    return image0    

#%%    

#%%   
#psf = GetDirtyBeam()
psf = psf_mat
#%%
#imsize = psf.shape[0] 
imsize = 512
dtbeam=fftshift(ifft2(ifftshift(psf),norm='ortho'))
#dtbeam=fftshift(ifft2(ifftshift(img.UVzeropad(psf,imsize)),norm='ortho'))
#dtbeam = fftshift(np.abs(ifft2(psf,[imsize,imsize],norm='ortho')))
#dtbeam = fftshift(np.real(ifft2(psf,[imsize,imsize],norm='ortho')))
fig1,(ax_b1)=plt.subplots(1,1)
ax_b1.imshow(dtbeam.real,cmap='jet')

mbx,mby=np.unravel_index(dtbeam.real.argmax(),(imsize,imsize))
kernel = dtbeam.real[mbx-3:mbx+3+1,mby-3:mby+3+1]

#%%

s = PointSource(imsize,3,(100,130),(200,200),(100,300))
noise = np.random.normal(0,1,(imsize,imsize))
map0 = s + noise

#kernel = kernel*100
map1=signal.fftconvolve(map0,dtbeam.real,mode='same')
#map1=signal.fftconvolve(map0,kernel,mode='same')
fig2, (ax_map,ax_dtmap) = plt.subplots(1,2)
ax_map.imshow(map0,cmap='jet')
ax_dtmap.imshow(map1,cmap='jet')

#%%

thresh = map1[300:400,300:400].std()
lpgain = 0.1
Niter = 600
clmodel = np.zeros_like(map1)
res = np.copy(map1)
#ros = int(np.floor(kernel.shape[0]/2))
#mx,my = 0, 0
for n in range(Niter):
    [mx,my] = np.unravel_index(np.fabs(res).argmax(),res.shape)
    clmodel[mx,my]+=lpgain*res[mx,my]
#    print 'mx=%d,my=%d'%(mx,my)
#    print 'max of res is %f'%(res[mx,my])
#    logger.debug('mx=%d,my=%d',mx,my)
#    logger.debug('max of res is %f', res[mx,my])
#    res[mx,my] = res[mx,my]-kernel.max()*res[mx,my]
    okxy,osxy=dv.effregion(kernel,res,mx,my)
    res[osxy[0]:osxy[1]+1,osxy[2]:osxy[3]+1] -= kernel[okxy[0]:okxy[1]+1,okxy[2]:okxy[3]+1]*lpgain*res[mx,my]
#    res[mx-ros:mx+ros+1,my-ros:my+ros+1] = 0.01*res[mx-ros:mx+ros+1,my-ros:my+ros+1]
#    logger.debug('After subtraction max of res is %f', res[mx,my])
    
    if (np.fabs(res).max()< thresh):
        break
#print 'After subtraction,max of res is %f'%(res[mx,my])    
#%%
## using gaussian beam to fit dirty beam
Imsize = 1*u.deg      
uvres = np.floor((1/Imsize.to(u.rad))).value
xg,yg = np.arange(-3,4),np.arange(-3,4)
xpg,ypg = np.meshgrid(xg,yg)

g_init = models.Gaussian2D(amplitude=1., x_mean=0, x_stddev=1.,
                           y_mean=0, y_stddev=1.,theta=0)
#fit_g = fitting.LevMarLSQFitter()
fitg = fitting.LevMarLSQFitter()
g2d = fitg(g_init, xpg, ypg,kernel)     
kg = g2d(xpg,ypg)                      
#plt.subplot(1,2,1)
#plt.imshow(g2d(xpg,ypg),cmap='jet')
#plt.subplot(1,2,2)
#plt.imshow(kernel,cmap='jet')
clmap = signal.fftconvolve(clmodel,kg,'same')
#%%
fig3, (ax_clmap,ax_res) = plt.subplots(1,2)
ax_clmap.imshow(map1,cmap='jet')
ax_res.imshow(clmap,cmap='jet')

#%%
def GetDirtyImage(imagemodel, dirtybeam):
    """
        Return a convolved dirtyimage    
    """

def HcleanTst(cleanbeam,dirtyimage):
    """
        Using a three-point-source to test Hogbom CLEAN      
    """

#%%
# delay reading test
import numpy as np
import muser_rdraw as mdr
import matplotlib.pyplot as plt
from _dlysys import dlysys
import rdin as ri
import cali
#%%
filename = r'C:\MUSER_DATA\CSRH_20141111-122131_98335183'
#fid = open(filename,'rb')
#offset = 321
#dly = mdr.dly_rd_m1(fid,offset)
#gpstime = mdr.gps_rd_m1(fid,offset)
#fid.close()

tobs = ri.getftr(filename)
trange = ['2014-11-11 12:22:45', '2014-11-11 12:22:55']
dout = ri.rdraw(filename,trange,tdec=5,array='m1', 
          nant=40,rfind=0,pol='LL')
dout0 = np.copy(dout['x'])
Nts = dout['time'].size
dlys = np.repeat(dlysys['2014-1111'][0:40].reshape(40,1),Nts,axis=1)
dout_fr = cali.frstop(dout,dlys)
#plt.plot(dout['x'][0:,3,0],'bo-')
#plt.plot(dout_fr['x'][0:,3,0],'rs-')
#%%
#obs_freq = 1.7025e+9 
##site_lat = 42.2 * np.pi / 180
##site_lon = 115.25 * np.pi / 180
#site_hr = (4*15) * np.pi / 180
#site_dec = 21.71*np.pi/180
#antflaglist = [13,26,39]
## Set imaging parameters
#Imdim = 512
##    uvdarr = np.zeros((Imdim,Imdim))
#Imsize = 1*u.deg      
#uvres = np.floor((1/Imsize.to(u.rad))).value
#rowt = 4*uvres
#
## need to flag antenna and make uvwbl conj symmetry
#uvwbl,flbllist = img.muvw('m1',antflaglist,obs_freq,site_hr,site_dec)
##plt.scatter(uvwbl['u'],uvwbl['v'],marker='.')
#
#
#gdpts = np.ceil(max(uvwbl['u'].max(),uvwbl['v'].max())/uvres)+1
#uvgspan = np.arange(-gdpts*uvres,(gdpts+1)*uvres,uvres)
#M = int(2*gdpts+1)
##    dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]
#Np = len(ma.compressed(uvwbl['u']))
#psf = np.array(np.zeros(Np),  dtype=ri.dtk)
#psf['x'] = np.ones(Np)
#psf['wgt'] = np.ones(Np)
#psf['u'] = ma.compressed(uvwbl['u'])
#psf['v']= ma.compressed(uvwbl['v'])
#psf['w']= ma.compressed(uvwbl['w'])
#
##%%
#psf_wt = img.mweight(psf,rowt)    
#
#deterx = np.ones(psf_wt['x'].shape[0]/2,dtype='bool')
#deteruv = np.ones((psf_wt['u'].shape[0]/2,2),dtype='bool')
#for i in range(psf_wt['x'].shape[0]/2):
#    deterx[i] = psf_wt['x'][i]==psf_wt['x'][i + psf_wt['x'].shape[0]/2]
#    deteruv[i,0] = psf_wt['u'][i]==-psf_wt['u'][i + psf_wt['u'].shape[0]/2]
#    deteruv[i,1] = psf_wt['v'][i]==-psf_wt['v'][i + psf_wt['v'].shape[0]/2]
#print deterx.all(), deteruv.all()
##%%
#kervspan = 3*uvres
#rms = uvres/2
#kerv = ('g',rms) # gaussian function is kernel
#psfg = img.mgrid(psf_wt,uvgspan,uvres,kerv,kervspan) 
#
##%%
#Npg=len(psfg)   
#deterxg = np.ones((Npg-1)/2+1,dtype='bool')
#deteruvg = np.ones(((Npg-1)/2+1,2),dtype='bool')
#for i in range((Npg-1)/2+1):
#    deterxg[i] = psfg['xg'][i]==np.conj(psfg['xg'][Npg-1-i])
#    deteruvg[i,0] = psfg['ug'][i]==-psfg['ug'][Npg-1-i]
#    deteruvg[i,1] = psfg['vg'][i]==-psfg['vg'][Npg-1-i]
#print deterxg.all(), deteruvg.all()
#
#plt.scatter(psfg['ug'],psfg['vg'],marker='.')
#indnz = np.where(np.logical_and(psfg['ug']==0,psfg['vg']==0))
##%%   
#psf_mat = img.mgrid2mat(psfg,M)
#
#determat = np.ones(psf_mat.shape,dtype='bool')
#for i in range(M):
#    for j in range(M):
#        determat[i,j]=(psf_mat[i,j]==np.conj(psf_mat[M-1-i,M-1-j]))
#
#beam = ifft2(ifftshift(psf_mat),norm='ortho')

#gmat = np.zeros(psf_mat.shape)
#gmat[M/2,M/2]=M*M
#gmat[M/2-10,M/2+11],gmat[M-1-(M/2-10),M-1-(M/2+11)]=M,M
    