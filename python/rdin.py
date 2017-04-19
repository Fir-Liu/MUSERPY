# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:48:23 2016

@author: Super
"""

import numpy as np
import muser_rdraw as mr
from astropy.time import Time
import misc as mi
import os
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
FRAME_NUM_FILE = 19200 # 19200 frames each file
FRAME_TICK = 0.003125 # 3.125ms for each frame
FRAMESIZE_MUSER_I = 100000
FRAMESIZE_MUSER_H = 204800
FILESIZE_MUSER_I = 1920000000 # filesize of a MUSER-I data file
FILESIZE_MUSER_H = 3932160000 # filesezie of a MUSER-H data file
HEADER = np.array(27*[85],dtype='u1')
#import datetime as dt
def get_file_tr(filename,array='m1'):
    """
    Usage:
        Read a raw muser data file, return the time range of the file.
    Returns:
        time_start, time_end
    """
    file_hd_addr = 0
    file_end_addr = FRAME_NUM_FILE - 1
    filesize = os.path.getsize(filename)
    fileframe = filesize//FRAMESIZE_MUSER_I if array=='m1' \
                else filesize//FRAMESIZE_MUSER_H
    with open(filename,'rb') as fid:
         if array == 'm1':
             fr_end_addr = file_end_addr if filesize == FILESIZE_MUSER_I \
                        else fileframe-1             
             t1 = mr.gps_rd_m1(fid,file_hd_addr)
             t2 = mr.gps_rd_m1(fid,fr_end_addr)       
         else:
             fr_end_addr = file_end_addr if filesize == FILESIZE_MUSER_H \
                        else fileframe-1             
             t1 = mr.gps_rd_m2(fid,file_hd_addr)
             t2 = mr.gps_rd_m2(fid,fr_end_addr)       
    return (t1,t2)                        
  
def get_db_tr(pathname,array='m1'):
    """
    Usage:
        get time-range of each datafile in database.
    Returns:
        [[fileA, (tA1, tA2)], [fileB, (tB1, tB2)]...]
    """    
    files = [os.path.join(pathname,n) for n in os.listdir(pathname) if os.path.isfile(os.path.join(pathname,n))]
    if array=='m1':
        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_I] #if filesize >1 frame it should be a datafile            
    else:
        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_H] #if filesize >1 frame it should be a datafile            
#    tranges = Time([get_file_tr(df,array) for df in datafiles])
    tfranges = [[df,get_file_tr(df,array)] for df in datafiles]
    return tfranges    

def dis_db_tr(pathname,array='m1'):
    """
    Usage:
        display the database's time-range.
    Returns:
        None
    """    
    tfranges = get_db_tr(pathname,array)
    tranges = Time([v[1] for v in tfranges])
    t_min = Time(np.min(tranges.jd),format='jd').isot
    t_max = Time(np.max(tranges.jd),format='jd').isot
    print('{}{}{}{}'.format('The database time range : ',t_min,' to: ',t_max))

def set_tr2db(pathname,array='m1'):    
    """
    Usage:
        set the default time-range as same as database time-range.
    Returns:
        time_start, time_end
    """    
    tfranges = get_db_tr(pathname,array)
    tranges = Time([v[1] for v in tfranges])
    t_min = Time(np.min(tranges.jd),format='jd').isot
    t_max = Time(np.max(tranges.jd),format='jd').isot    
    return [t_min, t_max]

def sel_file(pathname,trange,array='m1'):
    """
    Usage:
        determine which data-file should be processed in trange(time-range)
    Returns:
        filelist to be processed by rdraw() or others.
    """       
    tf = get_db_tr(pathname,array)
    tf1 = [[vt[0],Time(vt[1]).jd] for vt in tf]
    sortkf = lambda s: s[1][0]
    tfranges = sorted(tf1,key=sortkf)
    trange = Time(trange).jd
    trscope = np.array([v[1] for v in tfranges])
    filescope = [v[0] for v in tfranges]
    find1,find2 = -1,-1
#    find2 = -1
#    for ik,s in enumerate(trscope):
#        if trange[0]>=s[0] and trange[0]<=s[1]:
#            fileind.append(ik)
#            for im,v in enumerate(trscope[ik:]):
#                if trange[1]>=v[0] and trange[1]<=v[1]:
#                    fileind.append(im)
#                    break
#            break
    for ik,s in enumerate(trscope):
#        print(ik,s)
        if trange[0]>=s[0] and trange[0]<=s[1]:
            find1 = ik
#            print(ik,s)
#            print(fileind)
            break
        
#        if ik == 1: break
    for im,s in enumerate(trscope):
#        print(im,s)
        if trange[1]>=s[0] and trange[1]<=s[1]:
#            print('enter condition')
#            print(fileind)
            find2 = im
#            print(fileind)            
            break
#    print(fileind)
    if find1==-1 and find2==-1:
        print('time is out of datafile range, try to input again')
        return
    elif find1==-1:
#        inds = find2
        return filescope[find2]
    elif find2==-1:
        return filescope[find1]
    else:
        inds = list(range(find1,find2+1))
        sel_files = [filescope[p] for p in inds]
        return sel_files
#    return fileind            
#    return tfranges
    

def sorthelper(s):
    return s[1][0]
    
def rdraw(pathname,trange,ttick=0,array='m1', \
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
        ttick: type=float, time interval in seconds, ttick=0 means picking only 1
                time point data out,default = 0.
        array: type=str, array name in ['m1','m2'], default = 'm1'.
        nant:  type=int,number of antennas, default = 40.
        rfind: type=int, refered to which radio freq band is beginning in trange.
                for m1 rfind is in range(0,4); for m2 rfind is in range(0,33).
        pol: type=str,circular polarization in ['LL','RR'], default='L'.
    """              
    Nant,Nbl,Nfreq = nant,nant*(nant-1)/2,16
#    FRAME_TICK = 0.003125
    files = [os.path.join(pathname,n) for n in os.listdir(pathname) if os.path.isfile(os.path.join(pathname,n))]
    if(array == 'm2'):
        gps_rd = mr.gps_rd_m2
        rftag_rd = mr.rftag_rd_m2
        acor_rd = mr.acor_rd_m2
        xcor_rd = mr.xcor_rd_m2
        dly_rd = mr.dly_rd_m2
        print('Array is set to MUSER-II')
#        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_H] #if filesize >1 frame it should be a datafile    
    else:
        gps_rd = mr.gps_rd_m1
        rftag_rd = mr.rftag_rd_m1
        acor_rd = mr.acor_rd_m1
        xcor_rd = mr.xcor_rd_m1        
        dly_rd = mr.dly_rd_m1
        print('Array is set to MUSER-I')
#        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_I] #if filesize >1 frame it should be a datafile            
#   check time-range and select the right datafiles in time-range
#    for df in datafiles:
    tfranges = get_db_tr(pathname,array)

    tr = Time(trange)
#    logger.debug('tr is',tr)
    td = tr[1]-tr[0]
    if ttick == 0:
        Nts = 1
        fr_step = 0
    else:
        Nts = int(round(td.sec/ttick))
        fr_step = ttick/FRAME_TICK
    
    dout = {'x':np.zeros((Nbl,Nfreq,Nts),dtype='complex'), \
            'p':np.zeros((Nant,Nfreq,Nts)), \
            'GHz': np.zeros(Nts,dtype='S11'), \
            'time':np.zeros(Nts), \
            'dly': np.zeros((Nant,Nts))}
             
    bl2ord = mi.bl_tab(Nant)
    fid = open(filename,'rb')
    tini = Time(gps_rd(fid, 0))
    dt1 = tr[0]-tini
    fr_off = round(dt1.sec/FRAME_TICK)
    
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
#    
    print('Time range is:', trange)
    print('Number of Antenna:', nant)
    print('Poloarization:', pol)
    print('Rf band:', mi.ord2fr(rfind,array))
    print('{}{:3d}{:>10}'.format('Time bin interval:',ttick,'seconds'))
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
