# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 13:05:05 2016

@author: Super
"""

#%%        
#import rdin as ri
import misc as mi
import numpy.ma as ma
import img
import rdin as ri
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from numpy.fft import fft,ifft,fft2,ifft2,fftshift,ifftshift
from astropy import units as u
#from scipy import misc
#asc = misc.ascent()
#%%
#plt.clf()
#ft_asc = fft2(asc)
#asc_s = ifft2(fftshift(ft_asc))
#
#fig, (ax_asc,ax_ft_asc) = plt.subplots(1, 2)
#ax_asc.imshow(asc, cmap='gray')
#ax_ft_asc.imshow(asc_s.real, cmap='gray')
obs_freq = 1.7025e+9 
#site_lat = 42.2 * np.pi / 180
#site_lon = 115.25 * np.pi / 180
site_hr = (4*15) * np.pi / 180
site_dec = 21.71*np.pi/180
antflaglist = [12,13,25,26,38,39]
# Set imaging parameters
Imdim = 512
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
M = 2*gdpts+1
dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]
Np = len(ma.compressed(uvwbl['u']))
psf = np.array(np.zeros(Np),  dtype=dt)
psf['xk'] = np.ones(Np)
psf['wgtk'] = np.ones(Np)
psf['uk'] = ma.compressed(uvwbl['u'])
psf['vk']= ma.compressed(uvwbl['v'])
psf['wk']= ma.compressed(uvwbl['w'])

psf_wt = img.mweight(psf,rowt)    
kervspan = 3*uvres
rms = uvres/2
kerv = ('g',rms) # gaussian function is kernel
psfg = img.mgrid(psf_wt,uvgspan,uvres,kerv,kervspan)       

#%%
Nf = int(M)
Nz = Imdim**2 - Nf**2
Nzl,Nzr = Nz/2+1, Nz/2
psf_ex = np.hstack((np.zeros(Nzl),psfg['xg'],np.zeros(Nzr)))
psf_mat = psf_ex.reshape(Imdim,Imdim)
#    beam = ifft2(fftshift(psf_mat))
#    
#    beam = fftshift(np.real(beam))
beam = fftshift(np.abs(ifft2(psf_mat)))
plt.figure(4)
plt.imshow(beam)

#%%
psfgmat = psfg['xg'].reshape(255,255)
bmgmat = ifft2((psfgmat))
bmsgmat = ifft2(fftshift(psfgmat))
fig1,(ax_psf,ax_bm,ax_bms) = plt.subplots(1,3)
ax_psf.imshow(psfgmat,cmap='jet')
ax_bm.imshow(np.abs(bmgmat),cmap='jet')
ax_bms.imshow(np.abs(bmsgmat),cmap='jet')

res = bmsgmat-bmgmat
fig2,ax_res = plt.subplots(1,1)
ax_res.imshow(np.abs(bmsgmat-bmgmat),cmap='jet')
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
        image0[v[0]:v[0]+region+1,v[1]:v[1]+region+1] = 100 
    
    return image0    
imsize = 509
region = 8

s1 = PointSource(imsize,region,(100,200),(300,350),(200,150)) 
plt.imshow(s1,cmap='gray')

fs1 = ifftshift(np.pad(fftshift(fft2(s1)),(1,2),'constant',
                       constant_values=0))
s1z = ifft2(fs1)
plt.imshow(np.abs(s1z),cmap='gray')

fs1 = fftshift(fft2(s1))
plt.imshow(np.abs(fs1),cmap='gray')
z1 = np.zeros((511,1))
z2=np.zeros((1,512))
fs1z_ex=np.vstack((z2,np.hstack((z1,fftshift(fs1)))))
fs1z = fftshift(fs1z_ex)
s1z = ifft2(fs1z)
plt.imshow(np.abs(s1z),cmap='gray')
#%%
rpt= imsize-region
s = PointSource(imsize,region,(0,0),(0,rpt),(rpt,0),(rpt,rpt)) 

p = np.ceil(np.random.normal(0,3,(3,3)))
fp = fft2(p)
fps = fftshift(fp)

p1 = np.ceil(np.random.normal(0,4,(3,3)))
fp1 = fft2(p1)

z1 = np.zeros((3,1))
z2 = np.zeros((1,4))
fp1z_ex=np.vstack((z2,np.hstack((z1,fftshift(fp1)))))

pex = np.r_[np.zeros(4),np.ravel(p),np.zeros(3)].reshape(4,4)

fts = fft2(s)
fig1, (ax_s,ax_fts) = plt.subplots(1, 2)
ax_s.imshow(s, cmap='gray')
ax_fts.imshow(np.abs(fts), cmap='gray')   
#
fts1 = fftshift(fts)
s1 = ifft2(fts1)
fig2, (ax_s1,ax_fts1) = plt.subplots(1, 2)
ax_s1.imshow((np.abs(s1)), cmap='gray')
ax_fts1.imshow(np.abs(fts1), cmap='gray')   
#freqs = np.fft.fftfreq(9, d=1./9).reshape(3, 3)
#fs = ifftshift(freqs)
#%%
# test
kernel = signal.gaussian(16, 2)
ft_kernel = fft(kernel)
ft_kernel_zp = np.r_[ft_kernel,np.zeros(16)]
kernel_z = ifft(ft_kernel_zp)

ft_kernel_zp1=np.r_[np.zeros(8),fftshift(ft_kernel),np.zeros(8)]
kernel_z1 = ifft(ft_kernel_zp1)

fig1,(ax_ker,ax_ftker) = plt.subplots(2,1)
ax_ker.plot(kernel,'bo-')
ax_ftker.plot(ft_kernel.real,'b-')
ax_ftker.plot(ft_kernel.imag,'r-')

fig2,(ax_kerz,ax_ftkerz) = plt.subplots(2,1)
ax_kerz.plot(np.abs(kernel_z),'bo-')
ax_ftkerz.plot(ft_kernel_zp.real,'b-')
ax_ftkerz.plot(ft_kernel_zp.imag,'r-')                     

fig2,(ax_kerz1,ax_ftkerz1) = plt.subplots(2,1)
ax_kerz1.plot(np.abs(kernel_z1),'bo-')
ax_ftkerz1.plot(ft_kernel_zp1.real,'b-')
ax_ftkerz1.plot(ft_kernel_zp1.imag,'r-')                     

#%% test fftshift
kernel = signal.gaussian(15, 2)
ft_kernel = fft(kernel)
ft_kernel_s = fftshift(ft_kernel)
kernel_s = ifft(ft_kernel_s)

fig1,(ax_kernel,ax_ft_kernel) = plt.subplots(2,1)
ax_kernel.plot(kernel,'bo-')
ax_kernel.grid(True)
ax_ft_kernel.plot(ft_kernel.real,'bo-')
ax_ft_kernel.plot(ft_kernel.imag,'ro-')
ax_ft_kernel.grid(True)

fig2,(ax_kernel_s,ax_ft_kernel_s) = plt.subplots(2,1)
ax_kernel_s.plot(kernel_s.real,'bo-')
ax_kernel_s.plot(kernel_s.imag,'ro-')
ax_kernel_s.plot(abs(kernel_s),'k*-')
ax_kernel_s.grid(True)
ax_ft_kernel_s.plot(ft_kernel_s.real,'bo-')
ax_ft_kernel_s.plot(ft_kernel_s.imag,'ro-')
ax_ft_kernel_s.grid(True)

fig3,ax_kcp = plt.subplots(1,1)
ax_kcp.plot(kernel,'bo-')
#ax_kcp.plot(abs(kernel_s),'r*')
ax_kcp.plot(abs(kernel_s),'k*')
ax_kcp.grid(True)

#kernel_ex1 = np.r_[np.zeros(60),kernel]
#kernel_ex2 = np.r_[kernel,np.zeros(60)]
##x = np.arange(0,20)
#s = signal.bartlett(71)
#
#sc = signal.fftconvolve(s,kernel,'same')
#kernel_s = kernel_ex1 + kernel_ex2
#sc1 = signal.fftconvolve(s,kernel_s,'same')
#%%
fig1,(ax_kernel_s,ax_sc) = plt.subplots(2,1,sharex=True)
ax_kernel_s.plot(sc1,'r')
ax_kernel_s.grid(True)
#ax_kernel_s.set_xlim(0,20)
ax_sc.plot(sc)
ax_sc.grid(True)


#%%
#
#datafile = r'C:\MUSER_DATA\CSRH_20141111-122131_98335183'
#trange = ['2014-11-11 12:22:35', '2014-11-11 12:22:45']
#arrname,pol = 'm1','LL'
#Nant, rfind = 5, 0
#
#dout = ri.rdraw(datafile,trange,1,arrname,Nant,rfind,pol)
#
##%%
#import numpy as np
#from astropy import constants as const
#from astropy import units as u
#
#bl_min = 8*u.m
#bl_max = 3e3*u.m
#m1freq = np.array([0.4e9,2e9])*u.Hz
#m1wavlen = m1freq.to(u.m,equivalencies=u.spectral())
#m2freq = np.array([2e9,15e9])*u.Hz
#m2wavlen = m2freq.to(u.m,equivalencies=u.spectral())
#
#Imsizem1 = m1wavlen/bl_min*u.rad
#Imsizem2 = m2wavlen/bl_min*u.rad
#
#Imsize_set = 1*u.deg
#Impix_set = 512
#
#Bmsizem1 = m1wavlen/bl_max*u.rad
#Bmsizem2 = m2wavlen/bl_max*u.rad
#Impix = Imsize_set.to(u.rad)/(Bmsizem1/3)
#
#print 'M1 Imsize at 0.4 GHz (arcmin):', Imsizem1[0].to(u.arcmin).value
#print 'M1 Imsize at 2.0 GHz (arcmin):', Imsizem1[1].to(u.arcmin).value
#print 'M2 Imsize at 2.0 GHz (arcmin):', Imsizem2[0].to(u.arcmin).value
#print 'M2 Imsize at 15.0 GHz (arcmin):', Imsizem2[1].to(u.arcmin).value
#print 'M1 Bmsize at 0.4 GHz (arcsec):', Bmsizem1[0].to(u.arcsec).value
#print 'M1 Bmsize at 2.0 GHz (arcsec):', Bmsizem1[1].to(u.arcsec).value
#print 'M2 Bmsize at 2.0 GHz (arcsec):', Bmsizem2[0].to(u.arcsec).value
#print 'M2 Bmsize at 15.0 GHz (arcsec):', Bmsizem2[1].to(u.arcsec).value
##%%
import numpy as np
import scipy.signal as sig
import scipy.fftpack as sft
import matplotlib.pyplot as plt
#%%
#g = np.exp(-u**2/(2*s**2)
def gaussf(u,s):
    return np.exp(-u**2/(2*s**2))
            

x = np.arange(0,8)    
  
uk = np.array([1.2, 1.8, 2.06, 2.38, 3.5, 3.7, 4.9, 5.2])
data_uk = np.array(np.arange(1,9))
ukx = np.hstack((-uk[::-1],uk))
data_ukx = np.hstack((data_uk[::-1],data_uk))
win1 = gaussf(ukx,2)  
data_win_ukx = data_ukx*win1
gridwid = 1
kernelwidth = 4*gridwid
Kw = kernelwidth/gridwid
kernelsize = Kw + 1
ukgmax = np.ceil(ukx.max()/gridwid)
ukg = np.arange(-ukgmax,ukgmax+gridwid,gridwid)

#ukg_seg = np.concatenate([ukg[(len(ukg)-1)/2-Kw/2*gridwid:],np.arange(ukg[-1]+gridwid,ukg[-1]+(Kw/2+1)*gridwid)])
ukg_seg = np.arange(ukg[(len(ukg)-1)/2]-Kw/2*gridwid,ukg[-1]+(Kw/2+1)*gridwid,gridwid)


#%%

plt.plot(ukx,data_ukx,'bo',ukx,win1*data_ukx,'ro',ukx,win1,'go')
plt.grid(True)

data_toft = sft.fftshift(data_win_ukx)


#freq = sft.fftfreq(16,1e-6)
#
#x = np.arange(0,15)
#gbm = signal.gaussian(len(x), std=3)
#gbm1 = 1-gbm


#%%
#plt.figure(1)
#plt.plot(x,gbm,'b-',x,gbm1,'r-')
#plt.title(r"Gaussian gbm ($\sigma$=7)")
#plt.ylabel("Amplitude")
#plt.xlabel("Sample")

#%%
# 1-d u sampling and gridding simulation
u = np.array( [1.1, 1.3, 2.23,2.36,2.57,2.8, 3.13,3.26,3.78,4.08,4.32,
     4.57,4.8,5.25,5.68,6.25,7.83,9.41,10.6,13.8])
N = len(u)
r = np.ones(N)
win = lambda x: np.exp(-np.square(x)/40)
rwin = r*win(u)
rwin_hem = np.concatenate([rwin[::-1],rwin])
u_hem = np.concatenate([-u[::-1],u])
plt.stem(u_hem,rwin_hem)
plt.grid(True)
gdsize = 1.8
region = 4*gdsize
Lu = np.ceil(u.max()/gdsize)
gx = np.arange(-(Lu)*gdsize,(Lu+1)*gdsize,gdsize)
rg = np.ones(len(gx))
kspan = 4*gdsize
kernel = lambda x: np.exp(-x**2/8)
kx = np.linspace(-kspan/2,kspan/2,51)
plt.stem(gx,rg,'r-.'),plt.grid(True)
plt.plot(u_hem,rwin_hem,'k^')
plt.plot(kx,kernel(kx),'g-')
#%%
def grp(udata,xdata,start,stop,kspan,gdsize):
    grpdata = []
    kmax = (stop-start)/gdsize
    print kmax
    for k in np.arange(0,kmax+1):
        uarr = udata[np.where(np.logical_and(udata>start+k*gdsize-kspan/2,udata<start+k*gdsize+kspan/2))]
        rarr = xdata[np.where(np.logical_and(udata>start+k*gdsize-kspan/2,udata<start+k*gdsize+kspan/2))]
        grpdata.append(np.vstack((uarr,rarr)))
        print k,start+k*gdsize-kspan/2,start+k*gdsize+kspan/2
    return grpdata
#%%
grpdata = []
uarr = u_hem[np.where(np.logical_and(u_hem>-kspan/2,u_hem<kspan/2))]
rarr = rwin_hem[np.where(np.logical_and(u_hem>-kspan/2,u_hem<kspan/2))]
ur = np.vstack((uarr,rarr))
grpdata.append(ur)    
uarr1 = u_hem[np.where(np.logical_and(u_hem>gdsize-kspan/2,u_hem<gdsize+kspan/2))]
rarr1 = rwin_hem[np.where(np.logical_and(u_hem>gdsize-kspan/2,u_hem<gdsize+kspan/2))]
ur1 = np.vstack((uarr1,rarr1))
grpdata.append(ur1)

gdata = grp(u_hem,rwin_hem,0,Lu*gdsize,kspan,gdsize)

sv = np.zeros(len(gdata))
for k in range(len(gdata)):
    for j in range(len(gdata[k][0])):
        sv[k]=sv[k]+ kernel(gdata[k][0][j])*gdata[k][1][j]

#plt.figure(2)
#A = fft(gbm, 8) / (len(gbm)/2.0)
#A1= fft(gbm1, 8) / (len(gbm1)/2.0)
#freq = np.linspace(-0.5, 0.5, len(A))
#response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
#response1 = 20 * np.log10(np.abs(fftshift(A1 / abs(A1).max())))
#plt.plot(freq, response,'b-',freq, response1,'r-')
#plt.axis([-0.5, 0.5, -120, 0])
#plt.title(r"Frequency response of the Gaussian gbm ($\sigma$=7)")
#plt.ylabel("Normalized magnitude [dB]")
#plt.xlabel("Normalized frequency [cycles per sample]")
#
#%%
#2016/8/3
# flag list and numpy masked array
import numpy as np
import numpy.ma as ma
import misc as mc

bl2ord = mc.bl_tab(5)
uvdata = np.array([[2.1,2.2,2.3],[3.1,3.4,3.3],[4.2,5.6,5.6],[5.8,8.2,0],
                   [3.2,6.7,0],[2.8,4.3,0],[9.5,3.7,0],[11.1,10.3,0],
                    [13.4,12.7,0],[18.3,15.6,0]]).T
wgts = np.ones((1,10))
reg = 3.2
flaglist = [2]
muvdata = ma.masked_array(np.concatenate([uvdata,wgts]))
s1 = np.zeros((len(flaglist),5))
i = 0
for k in flaglist:
    s1[i,:] = bl2ord[k,:]
    i=i+1
s = s1[0,:]    
for j in range(s1.shape[0]-1):
    s = np.array(list(set(s)|set(s1[j+1,:])))
sind = np.array([k for k in s if k!=-1],dtype='u2')
for flind in sind:
    muvdata[:,flind]=ma.masked
    
cnt = 0    
for k in range(muvdata.shape[1]):
    if muvdata[0,k]>2:
       cnt = cnt + 1 
    else:       
        continue
cnt1 = 0        
for k in range(uvdata.shape[1]):
    if uvdata[0,k]>2:
       cnt1 = cnt1 + 1 
    else:       
        continue

for k in range(muvdata.shape[1]):
    for j in range(0,k)+range(k+1,muvdata.shape[1]):
        if (np.abs(muvdata[0,j]-muvdata[0,k])<=reg
            and np.abs(muvdata[1,j]-muvdata[1,k])<=reg):
           muvdata[3,k] = muvdata[3,k] + 1 
        else:       
            continue
#%%
#2016/8/2, 
#numpy array record or numpy structure array            
from scipy.spatial import distance as dis
s = np.array([[1,3,5,2,5],[2,1,-2,3,1]])  
s1 = s.T          
s2 = [tuple(u) for u in s1]
r = np.array([dis.euclidean(v,(0,0)) for v in s2],ndmin=2)
s3 = np.hstack([s1,r.T])
dt = [('uk',float),('vk',float),('xk',float),('wgtk',float)]
s4 = np.zeros((5,), dtype=dt)
s4['uk']=[u1 for u1 in s3[:,0]]
s4['vk']=[u2 for u2 in s3[:,1]]
s4['rk']=[u3 for u3 in s3[:,2]]
#s4 =np.array([values for values in s3],dtype=dt)

#%%    
#2016/8/3
#2016/8/9
# discrete convolution, weighting
import numpy as np
import misc as mc
import astropy.constants as con
from astropy import units as u
import matplotlib.pyplot as plt
import numpy.ma as ma
Nbl = 10
uk = np.around(10*np.random.randn(Nbl),1)  
vk = np.around(10*np.random.randn(Nbl),1)
wk = np.zeros(Nbl)
xk = np.ones(Nbl)
#%%
dt = [('uk',float),('vk',float),('wk',float),('xk',float),('wgtk',float)]  
xrec = np.zeros(2*Nbl,dtype=dt)
xrec['uk'] = np.hstack((uk,-uk))
xrec['vk'] = np.hstack((vk,-vk))
xrec['wk'] = np.hstack((wk,wk))
xrec['xk'] = np.hstack((xk,xk))
xrec['wgtk'] = np.ones(2*Nbl)

gdsize = 2.45
span = max(max(uk),max(vk))
gdpts = np.ceil(span/gdsize)+2



#%%
#2016/8/8, try improving 2-d discrete convolution performance
ker = lambda x,y : np.exp(-x**2/4-y**2/4)
kerspan = 1*gdsize
gdpts = 2
#gdpts = int(np.ceil(span/gdsize))+1
M = 2*gdpts+1
dt1 = [('ug',float),('vg',float),('xg',float)]
xconv = np.zeros(M**2,dtype=dt1)
uvgspan = np.arange(-gdpts*gdsize,(gdpts+1)*gdsize,gdsize)

#%%
def kerv(xseq, yseq):
    """
        return sequence: [f(xseq[0],yseq[0]),...]
    """
    kf = lambda x,y : np.exp(-x**2/4-y**2/4)    
    return map(kf,xseq,yseq)

def weight(psf,rowt): 
    # dtype is [('u',float),('v',float),('w',float),('x',float/complex),('wgt',float)]
    fdname = psf.dtype.names
    for idx,uvargs in enumerate(psf): # searching over half of psf is enough
        C1 = np.logical_and(np.abs(psf[fdname[0]]-uvargs[0])<rowt,psf[fdname[0]]!=uvargs[0])
        C2 = np.logical_and(np.abs(psf[fdname[1]]-uvargs[1])<rowt,psf[fdname[1]]!=uvargs[1])
        kind = np.where(np.logical_and(C1,C2))[0]
        psf[fdname[4]][idx] = len(kind) if len(kind)>0 else psf[fdname[4]][idx]
    return psf    
#%%
dt_wgt = [('uk',float),('vk',float),('xk',float),('wgtk',float)] 
xrec_wgt = np.zeros(len(xrec),dtype=dt_wgt)   
xrec_wgt['uk'],xrec_wgt['vk'] = xrec['uk'],xrec['vk']
xrec_wgt = weight(xrec_wgt,kerspan)    
#%%    
# Gridding with discrete convolution, stupid way 
i,j = 0,0 #  i:row index, j: column index
for ug in uvgspan:
    ind1 = np.where(np.abs(xrec['uk']-ug)<kerspan)[0]
    for vg in uvgspan[::-1]:
        ind2 = np.where(np.abs(xrec['vk']-vg)<kerspan)[0]
        kind = set(ind1)&set(ind2)
#        print kind,(ug,vg)
#        print (ug,vg)
        if len(kind)>0:
            xind = np.ravel_multi_index((i,j),(M,M))
#            print 'xind=',xind
#            xconv['xg'][xind]= -1            
            for kk in kind:
                xconv['xg'][xind] = xconv['xg'][xind] + \
                    xrec['xk'][kk]*ker(xrec['uk'][kk]-ug,xrec['vk'][kk]-vg)
        else:
            pass # do not use continue
        i = np.mod(i + 1,M)
    j = j + 1

xvs1 = np.copy(xconv['xg']).reshape((M,M))
#        return xrec
#%%
# Gridding with discrete convolution, clear way 
xconv1 = np.array([(ux,vx,0) for vx in uvgspan[::-1] for ux in uvgspan],dtype = dt1)

for idg,(ug,vg,_) in enumerate(xconv1):
    inda = np.where(np.logical_and(np.abs(xrec['vk']-vg)<kerspan,
                       np.abs(xrec['uk']-ug)<kerspan))[0]
    xconv1['xg'][idg] = xrec['xk'][inda].dot(kerv(xrec['uk'][inda]-ug, 
                        xrec['vk'][inda]-vg)) if len(inda)>0 else xconv1['xg'][idg]
#xvs2 = np.copy(xconv1['xg']).reshape((M,M))
#%%
# Gridding with discrete convolution, symmetry
xconv2 = np.array([(ux,vx,0) for vx in uvgspan[::-1] for ux in uvgspan],dtype = dt1)

for idg,(ug,vg,_) in enumerate(xconv2[0:(M**2-1)/2+1]):
    inda = np.where(np.logical_and(np.abs(xrec['vk']-vg)<kerspan,
                       np.abs(xrec['uk']-ug)<kerspan))[0]
    xconv2['xg'][idg] = xrec['xk'][inda].dot(kerv(xrec['uk'][inda]-ug, 
                        xrec['vk'][inda]-vg)) if len(inda)>0 else xconv2['xg'][idg]
    xconv2['xg'][M**2-1-idg] = xconv2['xg'][idg]
                        


#%%
span = np.arange(-1,2)
dtx = [('ux',float),('vx',float),('xx',float)]
#dtx = [('ux',float),('vx',float)]
uvx = np.array([(ux,vx) for vx in span[::-1] for ux in span ],dtype = dtx[0:2])
#uvx = (uvx.reshape((3,3)).T).flatten()
#ix,jx = 0, 0
#kx = 0
#for idx,(ux,vx,_) in enumerate(uvx):
#    uvx['xx'][idx] = uvx['xx'][idx] + idx
#    print idx
#%%
# plot rectangle
Nfig = gdpts+2
plt.figure(1)
plt.scatter(xrec['uk'],xrec['vk'],c='b',marker='.')
plt.scatter(0,0,c='r',marker='o')
plt.xticks([i*gdsize for i in range(-Nfig,Nfig+1)])
plt.yticks([i*gdsize for i in range(-Nfig,Nfig+1)])
plt.grid(True)
plt.xlim((-Nfig*gdsize,Nfig*gdsize))
plt.ylim((-Nfig*gdsize,Nfig*gdsize))

# rectangle plot
from matplotlib.patches import Rectangle
someX, someY = -gdpts*gdsize, -gdpts*gdsize
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((someX,someY), (M-1)*gdsize, (M-1)*gdsize, facecolor="none"))

#%%    
## loop estimation
M = np.array([2**i+1 for i in np.arange(11)])
n1 = M**2
n2= (M-1)*M/2+(M-1)/2+1
plt.figure(2)
plt.plot(M,n1,'r-')
plt.plot(M,n2,'b-')

    
#%%
#2016/8/5, try improving weighting performance
s1 = np.arange(1,5,1)    
s2 = s1-2.5
dt2 = [('ux',float),('vx',float),('wx',float),('x',float)]
xconv = np.zeros(len(s1),dtype=dt2)
xconv['ux']=s1
xconv['vx']=s2
#xconv-(xconv['ux'][0],xconv['vx'][0])

#np.array(xconv[0])-xconv['ux'][0]
#np.where(abs(xconv-(xconv['ux'][0],xconv['vx'][0])))
#np.where(xconv>np.array((xconv['ux'][0],xconv['vx'][0]),dtype=dt2))

for ux,vx,_,_ in xconv:
#    print ux,vx    
    kindu = set(np.where(np.abs(xconv['ux']-ux)<2)[0])
    kindzu = set(np.where(xconv['ux']==ux)[0])   
#    kindset = kindu-kindz if len(kindu)>0 else {}
    kindsetu = (kindu-kindzu) & kindu
    kindv = set(np.where(np.abs(xconv['vx']-vx)<2)[0])
    kindzv = set(np.where(xconv['vx']==vx)[0])   
#    kindset = kindu-kindz if len(kindu)>0 else {}
    kindsetv = (kindv-kindzv) & kindv
    kindset = kindsetu & kindsetv
    print kindset
#%%
# plot 3-d mesh
from matplotlib import cm
#from mpl_toolkits.mplot3d import Axes3D
ker = lambda x,y : np.exp(-x**2/4-y**2/4)
gdpts = np.ceil(span/gdsize)+2
xm = np.arange(-gdpts*gdsize,(gdpts+1)*gdsize,gdsize)
xx, yy = np.meshgrid(xm, xm, sparse=True)
z = ker(xx-5*gdsize,yy-5*gdsize)
fig2 = plt.figure()
#ax = fig2.gca(projection='3d')
norm = cm.colors.Normalize(vmax=abs(z).max(), vmin=-abs(z).max())
cmap = cm.PRGn
#h = plt.contourf(xm,xm,z)
cset1 = plt.contourf(xm, xm, z,cmap=cmap,norm=norm)
#surf = ax.plot_surface(xx, yy, z)

#%%
dtp1 = [('u',float)]
dtp2 = [('u',float),('x',float)]
p1 = ma.masked_array(np.arange(10),dtype=dtp1)  
p3 = np.array(np.zeros(len(ma.compressed(p1['u']))),dtype=dtp1)  
p1[-1] = ma.masked
#p2 = ma.masked_array(np.zeros(10),mask=p1.mask,dtype=dtp2)
#%%
s = np.array([3,5])
bl2ord = mi.bl_tab(4)
rclist = []
#bl2ord_flat = np.ravel(bl2ord)
for sk in list(s):
    indfl = np.where(bl2ord==sk)
    rowcol = zip(indfl[0],indfl[1])[0]
    rclist.append(rowcol)
    
#%%
#2016/08/22: shifted beam
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
from scipy.optimize import curve_fit
#%%
bm = np.array([0,0,1,2,1,0,0])
bm1 = shift(bm,3,cval=0)
plt.plot(bm,'bo-')
plt.plot(bm1,'r*-')
plt.grid(True)

#%% Gaussian Fitting
#Gaussian function
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def gauss_2d(x,y,a,cor,u1,u2,std1,std2):
    return a/(2*np.pi*std1*std2*(1-cor)**0.5)*np.exp(-(x-u1)**2)
# program

x = np.arange(-50,51.)
# sorry about these data!
y = gauss_function(x,1,0,10)
noise = 0.1*np.random.normal(0,1,len(x))
yn = y+noise
plt.plot(x,y,'b-') #Plot the curve, the gaussian is quite clear
plt.plot(x,yn,'ok') #Overplot the dots
plt.grid(True)

mean_e = np.mean(yn)
std_e = np.std(yn)
popt,pcov = curve_fit(gauss_function,x,yn,p0=[1,mean_e,std_e])
plt.plot(x,gauss_function(x,*popt),'r-')
popt1,pcov1 = curve_fit(gauss_function,x,yn)
plt.plot(x,gauss_function(x,*popt1),'g-')
#%%
# Try to fit data with 2-d gaussian
from scipy.stats import multivariate_normal
#%%
origin = 'lower'
x1, x2 = np.mgrid[-1:1:.01, -1:1:.01]
xv=np.dstack((x1,x2))
meanv= np.array([0,0])
covv = np.array([[2,0],[0,5]])

y = multivariate_normal(mean=meanv, cov=covv)
CS=plt.contourf(x1, x2, y.pdf(xv),10,
                  #[-1, -0.1, 0, 0.1],
                  #alpha=0.5,
                  cmap=plt.cm.bone,
                  origin=origin)

CS2 = plt.contour(CS, levels=CS.levels[::2],
                  colors='r',
                  origin=origin,
                  hold='on')
#%%
x = np.arange(-1,1,.01)
y = np.arange(-1,1,.01)
xv, yv = np.meshgrid(x, y)
pos = np.dstack((xv, yv))
rv = multivariate_normal([0.5, -0.2], [[2.0, 0.3], [0.3, 0.5]])
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.contourf(xv, yv, rv.pdf(pos))

#%%
x = np.arange(-5, 5, 0.1)
y = np.arange(-5, 5, 0.1)
xx, yy = np.meshgrid(x, y, sparse=True)
xg, yg = np.meshgrid(x, y)
z = np.sin(xx**2 + yy**2) / (xx**2 + yy**2)
z1 = np.sin(x**2 + y**2) / (x**2 + y**2)
h = plt.contourf(xg,yg,z)
#%%
#2016/08/23
#2-d gaussian function fitting- 512x512 matrix
import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from numpy import sin, cos, pi
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from scipy.ndimage.interpolation import shift
#%%
x = np.linspace(-1,1,512)
y = np.linspace(-1,1,512)
xg,yg = np.meshgrid(x,y)
pos = np.dstack((xg,yg))
meanv = [0,0]
sigmav = [[2.,0],[0,2.]]
rv = multivariate_normal.pdf(pos,meanv,sigmav)
bmsinc = sin(pi*xg)/(pi*xg)*sin(pi*yg)/(pi*yg)
#%%
plt.figure(1)
plt.subplot(1,2,1)
rvbeam = plt.imshow(rv)
plt.colorbar()
plt.subplot(1,2,2)
sincbeam = plt.imshow(bmsinc)
plt.colorbar()
#%%
x = np.linspace(-5,5,11)
y = np.linspace(-5,5,11)
xg,yg = np.meshgrid(x,y)
pos = np.dstack((xg,yg))
bmsinc = sin(pi*xg)/(pi*xg)*sin(pi*yg)/(pi*yg)

#%%
s= np.arange(0,20,2.25)
sk = 4.1
t0=np.floor((sk-s[0])/2.25)*2.25
t1=np.ceil((sk-s[0])/2.25)*2.25
#%%
x = np.linspace(-4,4,41)
bmsinc = np.sinc(x)


bmspan=np.where(bmsinc>np.max(bmsinc)/(2**0.5))
def gauss_function(x, a, x0, std):
    return a*np.exp(-(x-x0)**2/(2*std**2))
#x[bmspan]
popt1,pcov1 = curve_fit(gauss_function,x[bmspan],bmsinc[bmspan])

g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, x[bmspan], bmsinc[bmspan])

plt.plot(x,bmsinc,'r-')
plt.plot(x,bmsinc,'b*')
plt.grid(True)
plt.plot(x[bmspan],gauss_function(x[bmspan],*popt1),'g-')
plt.plot(x[bmspan],g(x[bmspan]),'ko')

#%%
import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from numpy import sin, cos, pi
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
from scipy.ndimage.interpolation import shift

def sinc2d(x,y):
    return np.sinc(x)*np.sinc(y)
#%%
x = np.linspace(-4,4,41)
y = np.linspace(-4,4,41)

xg,yg = np.meshgrid(x,y)
bmsinc = sinc2d(xg,yg)
#plt.imshow(bmsinc)
np.where(bmsinc > np.max(bmsinc)/(2**0.5))
pos=zip(np.where(bmsinc > np.max(bmsinc)/(2**0.5))[0],np.where(bmsinc > np.max(bmsinc)/(2**0.5))[1])
bmsinc[pos[0][0],pos[0][1]]
sinc2d(x[pos[0][1]],y[pos[0][0]])

#xypos = [(x[q],y[p]) for p,q in pos]
xpos = [x[qx] for _,qx in pos]
ypos = [y[px] for px,_ in pos]

xpg,ypg = np.meshgrid(xpos,ypos)

g_init = models.Gaussian2D(amplitude=1., x_mean=0, x_stddev=1.,
                           y_mean=0, y_stddev=1.,theta=0)
#fit_g = fitting.LevMarLSQFitter()
fitg = fitting.LinearLSQFitter()
g2d = fit_g(g_init, xpg, ypg,sinc2d(xpg,ypg))                           
plt.subplot(1,2,1)
plt.imshow(bmsinc)
plt.subplot(1,2,2)
plt.imshow(g2d(xg,yg))