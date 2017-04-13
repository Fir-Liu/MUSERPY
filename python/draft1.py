# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 10:44:21 2016

@author: Super
"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from numpy.fft import fft,ifft,fft2,ifft2,fftshift,ifftshift
from scipy import misc
picpath = r'..\pics'
def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])
#%%
#asc = misc.ascent()
cman = (mpimg.imread(picpath+'\\cameraman.tif'))
ftcman = fftshift(fft2(cman))
ftcmanex = ifftshift(np.pad(ftcman,(128,128),'constant',constant_values=0))
cmanex = ifft2(ftcmanex)

ftcman1 = fft2(cman)
ftcman1[0,0] = 0
cman1 = ifft2(ftcman1)
plt.figure(1)
plt.imshow(np.abs(cmanex),cmap='gray')
plt.figure(2)
plt.imshow(np.abs(cman1),cmap='gray')

#%%
def ongrid(data,dt,dim):
#    dt = data.dtype.names
    data_st = np.sort(data,order=[dt[0][0],dt[1][0]])
    datamat = np.zeros((dim,dim),dtype=dt[1][1])
    rcind_uv = [(i,j) for j in range(dim) for i in range(dim)[::-1]]
    for idx,(i,j) in enumerate(rcind_uv):
        datamat[i,j] = data_st[dt[1][0]][idx]                
    return datamat
                
s= (np.arange(9)+1).reshape(3,3)
#rcind = [(i,j) for i in range(3) for j in range(3)]
#uvind = [(ug,vg) for vg in range(1,-2,-1) for ug in range(-1,2)]         
dt = [('ug',float),('vg',float),('xg',complex)]
data = np.array([(0,0,7),(0,2.3,1+2j),(0,-2.3,1-2j),
                 (2.3,0,2+1j),(-2.3,0,2-1j),(-2.3,2.3,3+4j),
                 (2.3,-2.3,3-4j),(-2.3,-2.3,5-6j),(2.3,2.3,5+6j)],dtype=dt)      
               
data_st = np.sort(data,order=['ug','vg'])
rcind_uv = [(i,j) for j in range(3) for i in range(3)[::-1]]
datamat = np.zeros((3,3),dtype='complex')
for idx,(i,j) in enumerate(rcind_uv):
    print '%d-%d,(%f,%f)'%(i,j,data_st['ug'][idx],data_st['vg'][idx])
    datamat[i,j] = data_st['xg'][idx]

datamat1 = ongrid(data,dt,dim=3)

#%%
N=32
x = np.zeros(N)
x[N/2-2:N/2+2]=1

xs = fftshift(x)
fx = fftshift((fft(xs)))
fig1, (ax,afx)=plt.subplots(2,1)
ax.plot(xs,'ro')
ax.plot(xs,'r-')
afx.plot(fx.real,'bo')
afx.plot(fx.real,'b-')
afx.plot(fx.imag,'ro')
afx.plot(fx.imag,'r-')
