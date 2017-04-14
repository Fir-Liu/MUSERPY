# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:48:23 2016

@author: Super
"""

import numpy as np
import muser_rdraw as mr
from astropy.time import Time
import misc as mi
#%%
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
#logging.disable(logging.DEBUG)
#%%
dtk=[('u',float),('v',float),('w',float),('x',complex),('wgt',float)]
dtg =[('ug',float),('vg',float),('xg',complex)] 
dtuvw = [('u',float),('v',float),('w',float)]
#import datetime as dt
def getftr(filename,array='m1'):
    """
        Read a raw muser data file, return the time range of the file.
    """
    fr_head = 0
    fr_end = 19199
    fid = open(filename,'rb')    
    t1 = mr.gps_rd_m1(fid,fr_head) if array=='m1' else \
        mr.gps_rd_m2(fid,fr_head)
    t2 = mr.gps_rd_m1(fid,fr_end)if array=='m1' else \
        mr.gps_rd_m2(fid,fr_end)
    fid.close()
    return [t1,t2]


def rdraw(filename,trange,tdec=0,array='m1', \
          nant=40,rfind=0,pol='LL'):
    """ 
    Usage:
        Open one certain raw data file, and read data in, getting
        them stored in dictionary.
    Returns: dout, type=dict
        dout['x']: xcor data, type=complex array, shape=(Nbaseline,Nfreq=16,times)
        dout['p']: acor data, type=float array, shape=(Nant,Nfreq=16,times)
        dout['time']: gpstime data, type=float array,format=jd, shape=(times)
        dout['GHz']: radio frequency band, type=str array,like '0.8-1.2G'
    Input: 
        filename: type=str, datafile's path.
        trange: type=str list,like ['2014-11-11 12:22:35', '2014-11-11 12:22:45']
        tdec: type=float, time interval in seconds, tdec=0 means picking only 1
                time point data out,default = 0.
        array: type=str, array name in ['m1','m2'], default = 'm1'.
        nant:  type=int,number of antennas, default = 40.
        rfind: type=int, refered to which radio freq band is beginning in trange.
                for m1 rfind is in range(0,4); for m2 rfind is in range(0,33).
        pol: type=str,circular polarization in ['LL','RR'], default='L'.
    """              
    Nant,Nbl,Nfreq = nant,nant*(nant-1)/2,16
    tfr = 0.003125
    if(array == 'm2'):
        gps_rd = mr.gps_rd_m2
        rftag_rd = mr.rftag_rd_m2
        acor_rd = mr.acor_rd_m2
        xcor_rd = mr.xcor_rd_m2
        dly_rd = mr.dly_rd_m2
        print('Array is set to MUSER-II')
    else:
        gps_rd = mr.gps_rd_m1
        rftag_rd = mr.rftag_rd_m1
        acor_rd = mr.acor_rd_m1
        xcor_rd = mr.xcor_rd_m1        
        dly_rd = mr.dly_rd_m1
        print('Array is set to MUSER-I')    
    tr = Time(trange)
#    logger.debug('tr is',tr)
    td = tr[1]-tr[0]
    if tdec == 0:
        Nts = 1
        fr_step = 0
    else:
        Nts = int(round(td.sec/tdec))
        fr_step = tdec/tfr
    
    dout = {'x':np.zeros((Nbl,Nfreq,Nts),dtype='complex'), \
            'p':np.zeros((Nant,Nfreq,Nts)), \
            'GHz': np.zeros(Nts,dtype='S11'), \
            'time':np.zeros(Nts), \
            'dly': np.zeros((Nant,Nts))}
             
    bl2ord = mi.bl_tab(Nant)
    fid = open(filename,'rb')
    tini = Time(gps_rd(fid, 0))
    dt1 = tr[0]-tini
    fr_off = round(dt1.sec/tfr)
    
    freq_list = np.zeros(9,dtype='u1')
    for ss in range(0,9):
        freq_list[ss] = mi.fr2ord(rftag_rd(fid, fr_off+ss))
    
    rfind1 = np.array([rfind,rfind],dtype='u1')
    pos = -1
    for k in range(0,len(freq_list)-1): 
        if (np.equal(rfind1,freq_list[k:k+2]).all()):
            pos = k
            break
        else:
            continue
    #tst = Time(mr.gps_rd_m1(fid, fr_off))
    if pos == -1:
        print('the observing mode should be no jumping, rfind \
                is neglected')    
    else:        
        if pol == 'LL':
            fr_off = fr_off+pos
        else:
            fr_off = fr_off+pos+1

    logger.debug('Frame starts from: %d',fr_off)
    logger.debug('Frame steps: %d',fr_step)        
    
    print('Time range is:', trange)
    print('Number of Antenna:', nant)
    print('Poloarization:', pol)
    print('Rf band:', mi.ord2fr(rfind,array))
    print('{}{:3d}{:>10}'.format('Time bin interval:',tdec,'seconds'))
    print('Time bins:', Nts)
    print('Waiting....................')
    for nn in range(Nts):
#        dout['time'][nn] = Time(mr.gps_rd_m1(fid, fr_off+nn*fr_step)).jd \
#                            if (array=='m1') else Time(mr.gps_rd_m2(fid, fr_off+nn*fr_step)).jd
#        dout['GHz'][nn] = mr.rftag_rd_m1(fid,fr_off+nn*fr_step) \
#                            if (array=='m1') else mr.rftag_rd_m2(fid,fr_off+nn*fr_step)
        dout['time'][nn] = Time(gps_rd(fid, fr_off+nn*fr_step)).jd 
        dout['GHz'][nn] = rftag_rd(fid,fr_off+nn*fr_step)
        dout['dly'][:,nn] = dly_rd(fid, fr_off+nn*fr_step)
        for aa in range(Nant):
#            dout['p'][aa,:,nn] = mr.acor_rd_m1(fid, fr_off+nn*fr_step,aa) \
#                                if(array=='m1') else mr.acor_rd_m2(fid, fr_off+nn*fr_step,aa)
            dout['p'][aa,:,nn] = acor_rd(fid, fr_off+nn*fr_step,aa)
            
            for kk in range(Nfreq):
#                dout['x'][bl2ord[aa,aa+1:],kk,nn]=mr.xcor_rd_m1(fid, fr_off+nn*fr_step, kk, aa,range(aa+1,Nant)) \
#                                                if(array=='m1') else mr.xcor_rd_m2(fid, fr_off+nn*fr_step, kk, aa,range(aa+1,Nant))
#                print('nn=',nn)
#                print('aa=',aa)
#                print('kk=',kk)
#                print('fr_off=',fr_off)
#                print('fr_step=',fr_step)
                dout['x'][bl2ord[aa,aa+1:],kk,nn]=xcor_rd(fid, fr_off+nn*fr_step, kk, aa,list(range(aa+1,Nant))) 
    
    fid.close()
    print('Reading data is done.')             
    return dout

#%% Test


#%%
