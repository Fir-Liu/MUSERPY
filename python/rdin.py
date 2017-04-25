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
        (time_start, time_end)
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
        get time-range of each datafile in database in ascending order
    Returns:
        [[fileA, (tA1, tA2)], [fileB, (tB1, tB2)]...]
    """    
    files = [os.path.join(pathname,n) for n in os.listdir(pathname) if os.path.isfile(os.path.join(pathname,n))]
    if array=='m1':
        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_I] #if filesize >1 frame it should be a datafile            
    else:
        datafiles = [s for s in files if os.path.getsize(s)>FRAMESIZE_MUSER_H] #if filesize >1 frame it should be a datafile            
#    tranges = Time([get_file_tr(df,array) for df in datafiles])
    tf = [[df,get_file_tr(df,array)] for df in datafiles]
    sortkf = lambda s: s[1][0]
    tfranges = sorted(tf,key=sortkf)
    return tfranges    

def dis_db_tr(pathname,array='m1'):
    """
    Usage:
        display the database's time-range.
    Returns:
        None
    """    
    t_min,t_max = set_tr2db(pathname,array)
    print('{}{}{}{}'.format('The database time range : \n',t_min,' to: ',t_max))

def set_tr2db(pathname,array='m1'):    
    """
    Usage:
        set the default time-range as same as database time-range.
    Returns:
        time_start, time_end
    """    
    tfranges = get_db_tr(pathname,array)
#    sortkf = lambda s: Time(s[1][0]).jd
#    tfranges = sorted(tf,key=sortkf)
    t_min,t_max = tfranges[0][1][0], tfranges[-1][1][1]
    return [t_min, t_max]

def set_trange(trange,pathname,array='m1'):
    """
    Usage:
        set the new time range based on the time-range of database files.
    Returns:
        time_start, time_end
    """     
    dbtf = get_db_tr(pathname,array)
    dbtfjd = [[vt[0],Time(vt[1]).jd] for vt in dbtf]
    trscope = np.array([v[1] for v in dbtfjd])   
    tr = Time(trange).jd
    for ik, s in enumerate(trscope):
        if ik == 0:
            tr[0] = s[0] if tr[0]<s[0] else tr[0] 
        else:
            if tr[0]>trscope[ik-1][1] and tr[0]<s[0]:
                tr[0] = s[0]
    for ik, s in enumerate(trscope):
        if ik == len(trscope)-1:
            tr[1] = s[1] if tr[1]>s[1] else tr[1] 
        else:
            if tr[1]<trscope[ik+1][0] and tr[0]>s[1]:
                tr[1] = s[1]                
    return tr

def sel_file(pathname,trange,array='m1'):
    """
    Usage:
        determine which data-file should be processed in trange(time-range)
    Returns:
        [fileA, fileB,...]
    """       
    trange = set_trange(trange,pathname,array) # respecify the input trange
    dbtf = get_db_tr(pathname,array)
    dbtfjd = [[vt[0],Time(vt[1]).jd] for vt in dbtf]

    trscope = np.array([v[1] for v in dbtfjd])
    filescope = [v[0] for v in dbtfjd]
    find1,find2 = -1,-1 # find1 start_index, find2: end_index
    for ik,s in enumerate(trscope):
        if trange[0]>=s[0] and trange[0]<=s[1]:
            find1 = ik
            break
    for im,s in enumerate(trscope):
        if trange[1]>=s[0] and trange[1]<=s[1]:
            find2 = im
            break
    if find1==-1 and find2==-1:
        print('time is out of datafile range, try to input again')
        return
    elif find1==-1:
        return filescope[find2]
    elif find2==-1:
        return filescope[find1]
    else:
        inds = list(range(find1,find2+1))
        sel_files = [filescope[p] for p in inds]
        return sel_files   

def count_file_addr(filename,array='m1'):
    filesize = os.path.getsize(filename)
    framenum = filesize//FRAMESIZE_MUSER_I if array=='m1' \
                else FRAMESIZE_MUSER_H
    return framenum

def get_files_addr(sel_files,trange,array='m1'):
    """
    Usage:
        determine start_frame_offset and end_frame_offset of
        given files and trange.
    Returns:
        [[fileA,(addrA1,addrA2)],[fileB,(addrB1,addrB2)]...]
    """     
    tobj = Time(trange)
    gps_rd = mr.gps_rd_m1 if array=='m1' else mr.gps_rd_m2   
    if isinstance(sel_files,list):
        addrlist = len(sel_files)*['']
        with open(sel_files[0],'rb') as fid1:
            t1=Time(gps_rd(fid1,0))
            dt1 = tobj[0]-t1
#            print(dt1.sec,FRAME_TICK)
            addrlist[0] =(int(mi.qround(dt1.sec,FRAME_TICK)/FRAME_TICK), \
                    count_file_addr(sel_files[0],array)-1)
#            print(addrlist[0])
        with open(sel_files[-1],'rb') as fid2:
            t2=Time(gps_rd(fid2,0))
            dt2 = tobj[1]-t2
            addrlist[-1] =(0,int(mi.qround(dt2.sec,FRAME_TICK)/FRAME_TICK))
        for ik, fs in enumerate(sel_files):
            if ik == 0 or ik == len(sel_files)-1:
                continue
            else:
                addrlist[ik] = (0,count_file_addr(sel_files[ik],array)-1)                       
        return addrlist
    else:
        with open(sel_files,'rb') as fid:
            t0 = Time(gps_rd(fid,0))
            dt1 = tobj[0]-t0
            dt2 = tobj[1]-t0
            fr1 = round(dt1.sec/FRAME_TICK)
            fr2 = dt2.sec//FRAME_TICK-1
        return fr1,fr2
                      
    
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
    Nant,Nbl,Nfreq = nant,int(nant*(nant-1)/2),16
    gps_rd = mr.gps_rd_m1 if array == 'm1' else mr.gps_rd_m2
    rftag_rd = mr.rftag_rd_m1 if array == 'm1' else mr.rftag_rd_m2
    acor_rd = mr.acor_rd_m1 if array == 'm1' else mr.acor_rd_m2
    xcor_rd = mr.xcor_rd_m1 if array == 'm1' else mr.xcor_rd_m2       
    dly_rd = mr.dly_rd_m1 if array == 'm1' else mr.dly_rd_m2
    arrayname = 'MUSER-I' if array =='m1' else 'MUSER-H'
    print('Array is set to '+arrayname)
#   check time-range and select the right datafiles in time-range
    bl2ord = mi.bl_tab(Nant)
    dfiles = sel_file(pathname,trange,array)
    addr_dfiles = get_files_addr(dfiles,trange)
    fr_step = int(ttick/FRAME_TICK) 
    if isinstance(dfiles,list):
        Ntsl = len(dfiles)*[0]
        for ik,a in enumerate(addr_dfiles):
            Ntsl[ik] =  1 if (a[1]-a[0]+1)//fr_step == 0 \
                        else (a[1]-a[0]+1)//fr_step
    else:
        Nts = [1] if (addr_dfiles[1]-addr_dfiles[0]+1)//fr_step == 0 \
                    else [(addr_dfiles[1]-addr_dfiles[0]+1)//fr_step]
    Nts = sum(Ntsl)
    print( fr_step,Nts,Ntsl )
    print(addr_dfiles); return
    dout = {'x':np.zeros((Nbl,Nfreq,Nts),dtype='complex'), \
            'p':np.zeros((Nant,Nfreq,Nts)), \
            'GHz': np.zeros(Nts,dtype='S11'), \
            'time':np.zeros(Nts), \
            'GPS':np.zeros(Nts,dtype='S26'), \
            'dly': np.zeros((Nant,Nts))}

    print('Time range is:', trange)
    print('Number of Antenna:', nant)
    print('Poloarization:', pol)
    print('Rf band:', mi.ord2fr(rfind,array))
    print('{}{:3d}{:>10}'.format('Time bin interval:',ttick,'seconds'))
    print('Time bins:', Nts)
    print('Waiting....................')
    
    for i,fn in enumerate(dfiles):
        with open(fn,'rb') as fid:
           for j in range(Ntsl[i]):
#               print(i,j)
               dout['GHz'][j + i*Ntsl[i-1]] = rftag_rd(fid, addr_dfiles[i][0] + j*fr_step)
               dout['GPS'][j + i*Ntsl[i-1]] = gps_rd(fid, addr_dfiles[i][0] + j*fr_step)
               for aa in range(Nant):
                   dout['p'][aa,:,j + i*Ntsl[i-1]] = acor_rd(fid, addr_dfiles[i][0] + j*fr_step,aa)
                   for kk in range(Nfreq):
                       dout['x'][bl2ord[aa,aa+1:],kk,j + i*Ntsl[i-1]]=xcor_rd(fid, addr_dfiles[i][0] + j*fr_step,kk, aa,list(range(aa+1,Nant)))
    print('Reading data is done.')  
    return dout
     
#    bl2ord = mi.bl_tab(Nant)
#    fid = open(filename,'rb')
#    tini = Time(gps_rd(fid, 0))
#    dt1 = tr[0]-tini
#    fr_off = round(dt1.sec/FRAME_TICK)
#    
#    freq_list = np.zeros(9,dtype='u1')
#    for ss in range(0,9):
#        freq_list[ss] = mi.fr2ord(rftag_rd(fid, fr_off+ss))
#    
#    rfind1 = np.array([rfind,rfind],dtype='u1')
#    pos = -1
#    for k in range(0,len(freq_list)-1): 
#        if (np.equal(rfind1,freq_list[k:k+2]).all()):
#            pos = k
#            break
#        else:
#            continue
#    #tst = Time(mr.gps_rd_m1(fid, fr_off))
#    if pos == -1:
#        print('the observing mode should be no jumping, rfind \
#                is neglected')    
#    else:        
#        if pol == 'LL':
#            fr_off = fr_off+pos
#        else:
#            fr_off = fr_off+pos+1

#    logger.debug('Frame starts from: %d',fr_off)
#    logger.debug('Frame steps: %d',fr_step)        
##    
#    print('Time range is:', trange)
#    print('Number of Antenna:', nant)
#    print('Poloarization:', pol)
#    print('Rf band:', mi.ord2fr(rfind,array))
#    print('{}{:3d}{:>10}'.format('Time bin interval:',ttick,'seconds'))
#    print('Time bins:', Nts)
#    print('Waiting....................')
#    for nn in range(Nts):
##        dout['time'][nn] = Time(mr.gps_rd_m1(fid, fr_off+nn*fr_step)).jd \
##                            if (array=='m1') else Time(mr.gps_rd_m2(fid, fr_off+nn*fr_step)).jd
##        dout['GHz'][nn] = mr.rftag_rd_m1(fid,fr_off+nn*fr_step) \
##                            if (array=='m1') else mr.rftag_rd_m2(fid,fr_off+nn*fr_step)
#        dout['time'][nn] = Time(gps_rd(fid, fr_off+nn*fr_step)).jd 
#        dout['GHz'][nn] = rftag_rd(fid,fr_off+nn*fr_step)
#        dout['dly'][:,nn] = dly_rd(fid, fr_off+nn*fr_step)
#        for aa in range(Nant):
##            dout['p'][aa,:,nn] = mr.acor_rd_m1(fid, fr_off+nn*fr_step,aa) \
##                                if(array=='m1') else mr.acor_rd_m2(fid, fr_off+nn*fr_step,aa)
#            dout['p'][aa,:,nn] = acor_rd(fid, fr_off+nn*fr_step,aa)
#            
#            for kk in range(Nfreq):
##                dout['x'][bl2ord[aa,aa+1:],kk,nn]=mr.xcor_rd_m1(fid, fr_off+nn*fr_step, kk, aa,range(aa+1,Nant)) \
##                                                if(array=='m1') else mr.xcor_rd_m2(fid, fr_off+nn*fr_step, kk, aa,range(aa+1,Nant))
##                print('nn=',nn)
##                print('aa=',aa)
##                print('kk=',kk)
##                print('fr_off=',fr_off)
##                print('fr_step=',fr_step)
#                dout['x'][bl2ord[aa,aa+1:],kk,nn]=xcor_rd(fid, fr_off+nn*fr_step, kk, aa,list(range(aa+1,Nant))) 
#    
#    fid.close()
#    print('Reading data is done.')             
#    return dout

#%% Test


#%%
