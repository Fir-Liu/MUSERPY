import numpy as np

def xcor_rd_m1( fid, frame_offset, chan, antref, ant_list ):
    '''
    Output:
        Return the complex xcorrelation data of baseline: antref--ant_list[k]
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
        chan: 0~15, channel number in 400 MHz
        antref: referenced antenna. 0~43
        antlist: [ant1, ant2, ant3...] 0~43
    '''
    zone_size = 12128
    block_size = 12
    slot_size = 2
    
    xcor_I_out = np.zeros((len(ant_list),1))
    xcor_Q_out = np.zeros((len(ant_list),1))
    dt = 'u1' #uint8
    xtable = np.zeros((44,44))
    xtable[0,:] = np.arange(44) #0:43;

    for k in range(43): #0:42
        xtable[k+1,:] = xtable[k,:] + 43-(k+1)

    zone_ind = np.ceil((chan+1)/2.)
    slot_ind = np.array([1, 3, 5]) if (np.mod((chan+1), 2) == 0) \
                                    else np.array([2,4,6])

    for n in range(len(ant_list)):             #1:length(ant_list)
        
#        print('n=',n)
        block_ind = xtable[antref,ant_list[n]];
        xcor_addr = frame_offset*100000 + 2944 + (zone_ind - 1)*zone_size \
                    +(block_ind - 1)*block_size + (slot_ind - 1)*slot_size    
                            
        xcor_addr = xcor_addr.astype(int)         
#        print('xcor_addr=',xcor_addr)
        
        fid.seek(0,0); fid.seek(xcor_addr[0],1)                                
        xcor_Q1 = np.fromfile(fid, dtype=dt,count=1)                
        xcor_I1 = np.fromfile(fid, dtype=dt,count=1)
        fid.seek(0,0); fid.seek(xcor_addr[1],1)
        xcor_Q2 = np.fromfile(fid, dtype=dt,count=1)                
        xcor_I2 = np.fromfile(fid, dtype=dt,count=1)        
        fid.seek(0,0); fid.seek(xcor_addr[2],1)
        xcor_Q3 = np.fromfile(fid, dtype=dt,count=1)                
        xcor_I3 = np.fromfile(fid, dtype=dt,count=1)        

        xcor_I_out[n] = xcor_I3*(256**2) + xcor_I2*256 + xcor_I1
        xcor_Q_out[n] = xcor_Q3*(256**2) + xcor_Q2*256 + xcor_Q1        
        xcor_I_out[n] = xcor_I_out[n]-2**24 if(xcor_I_out[n]>=2**23) else xcor_I_out[n]
        xcor_Q_out[n] = xcor_Q_out[n]-2**24 if(xcor_Q_out[n]>=2**23) else xcor_Q_out[n]        
        
    fid.seek(0, 0)        
    return list(map(np.complex, xcor_I_out, xcor_Q_out))

def acor_rd_m1(fid, frame_offset, ant_ind):
    '''
    Output:
        Return the autocorrelation data of 16 channels of one antenna in Muser-I
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
        ant_ind: 0~43
    '''
    P = np.zeros(16)
    dt = np.dtype('u4') #uint32
    tablist = np.hstack((np.reshape(np.arange(1,45),(44,1)),\
            np.reshape(np.repeat(np.reshape(np.arange(1,12),(1,11)),4,axis=0),\
            (44,1),order='F'),np.reshape(np.repeat(np.reshape(np.arange(1,5),\
            (4,1)),11,axis=1),(44,1),order='F')))
    
    frame_addr = frame_offset*100000 
    auto_data_addr = 2944 + 11392 
    zone_size = 12128
    block_size = 64
    slot_size = 8
    
    for kk in np.arange(8)+1:
        ant_auto_addr = frame_addr + auto_data_addr + (kk-1)*zone_size \
        +(tablist[ant_ind, 1] - 1)*block_size + (tablist[ant_ind, 2] - 1) \
        *slot_size
        ant_auto_addr = ant_auto_addr.astype(int)
        fid.seek(0,0)
        fid.seek(ant_auto_addr,1)
        P[2*kk-2]=np.fromfile(fid,dtype=dt,count=1)
        P[2*kk-1]=np.fromfile(fid,dtype=dt,count=1)
#        print ant_auto_addr
    fid.seek(0,0)
    return P

def rftag_rd_m1(fid, frame_offset):
    '''
    Output:
        Return the rf frequency label(string) of one frame in Muser-I datafile
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
    '''
    f_dict = {0:'0.4G-0.8G',85:'0.8G-1.2G',170:'1.2G-1.6G',255:'1.6G-2.0G'}
    dt = np.dtype('u1') #uint8
    rftag_addr = frame_offset*100000 + 208
    rftag_addr1 = frame_offset*100000 + 99264
    fid.seek(rftag_addr, 0)
    tagA = np.zeros(2)
    tagA[0] = np.fromfile(fid, dtype=dt,count=1)
    fid.seek(rftag_addr1,0)       
    tagA[1] = np.fromfile(fid, dtype=dt,count=1)
    fid.seek(0,0)
    return f_dict[tagA[1]]

def gps_rd_m1(fid, frame_offset):
    '''
    Output:
        Return the GPStime(string:isot format) of one frame in Muser-I file
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
    '''
    dt = np.dtype('u1') #uint8
    gps_addr = frame_offset*100000 + 32
    gps_info = np.zeros(8)    
    
    fid.seek(gps_addr,0)
    gps_info = np.fromfile(fid, dtype=dt, count=8)
    fid.seek(0,0)
    gps_st = {'yr':0,'mon':0,'day':0,'hr':0,'min':0,'sec':0,\
            'msec':0,'usec':0,'nsec':0}    
    s=np.array(list(map(np.binary_repr,gps_info,[8,8,8,8,8,8,8,8])))
    
#    gps_str = s[::-1].tostring()
    gps_str = ''.join(s[::-1])

    str1 = [gps_str[0:12],gps_str[12:16],gps_str[16:21],gps_str[21:26],gps_str[26:32],\
            gps_str[32:38],gps_str[38:48],gps_str[48:58],gps_str[58:64]]
    
    gps_st['yr'], gps_st['mon'],gps_st['day'],gps_st['hr'],gps_st['min'], \
    gps_st['sec'],gps_st['msec'],gps_st['usec'],gps_st['nsec']  \
    = list(map(int,str1,[2,2,2,2,2,2,2,2,2]))
    
    gps_st['yr']=gps_st['yr']+2000    

    gps_isot = '{yr}-{mon:02}-{day:02}T{hr:02}:{min:02}:{sec:02}.{msec:03}{usec:03}' \
            .format(**gps_st)
    return gps_isot

def gps_rd_m2(fid, frame_offset):
    '''
    Output:
        Return the GPStime(string:isot format) of one frame in Muser-II file
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
    '''
    dt = np.dtype('u1') #uint8
    gps_addr = frame_offset*204800 + 32
    gps_info = np.zeros(8)    
    
    fid.seek(gps_addr,0)
    gps_info = np.fromfile(fid, dtype=dt, count=8)
    fid.seek(0,0)
    gps_st = {'yr':0,'mon':0,'day':0,'hr':0,'min':0,'sec':0,\
            'msec':0,'usec':0,'nsec':0}    
    s=np.array(list(map(np.binary_repr,gps_info,[8,8,8,8,8,8,8,8])))
    
#    gps_str = s[::-1].tostring()
    gps_str = ''.join(s[::-1])

    str1 = [gps_str[0:12],gps_str[12:16],gps_str[16:21],gps_str[21:26],gps_str[26:32],\
            gps_str[32:38],gps_str[38:48],gps_str[48:58],gps_str[58:64]]
    
    gps_st['yr'], gps_st['mon'],gps_st['day'],gps_st['hr'],gps_st['min'], \
    gps_st['sec'],gps_st['msec'],gps_st['usec'],gps_st['nsec']  \
    = list(map(int,str1,[2,2,2,2,2,2,2,2,2]))
    
    gps_st['yr']=gps_st['yr']+2000    

    gps_isot = '{yr}-{mon:02}-{day:02}T{hr:02}:{min:02}:{sec:02}.{msec:03}{usec:03}' \
            .format(**gps_st)
    return gps_isot


def rftag_rd_m2(fid, frame_offset):
    '''
    Output:
        Return the rf frequency label(string) of one frame in Muser-II datafile
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
    '''
    f_dict = {0:'0.4G-0.8G',1:'0.8G-1.2G',2:'1.2G-1.6G',3:'1.6G-2.0G',\
            4:'0.4G-0.8G',5:'0.8G-1.2G',6:'1.2G-1.6G',7:'1.6G-2.0G',\
            8:'0.4G-0.8G',9:'0.8G-1.2G',10:'1.2G-1.6G',11:'1.6G-2.0G', \
            12:'0.4G-0.8G',13:'0.8G-1.2G',14:'1.2G-1.6G',15:'1.6G-2.0G', \
            16:'0.4G-0.8G',17:'0.8G-1.2G',18:'1.2G-1.6G',19:'1.6G-2.0G', \
            20:'0.4G-0.8G',21:'0.8G-1.2G',22:'1.2G-1.6G',23:'1.6G-2.0G', \
            24:'0.4G-0.8G',25:'0.8G-1.2G',26:'1.2G-1.6G',27:'1.6G-2.0G', \
            28:'0.4G-0.8G',29:'0.8G-1.2G',30:'1.2G-1.6G',31:'1.6G-2.0G', \
            32:'0.4G-0.8G'}
    dt = np.dtype('u1') #uint8
    rftag_addr = frame_offset*204800 + [208, 204664]
    fid.seek(rftag_addr[0], 0)
    tagA = np.zeros(2)
    tagA[0] = np.fromfile(fid, dtype=dt,count=1)
    fid.seek(rftag_addr[1],0)       
    tagA[1] = np.fromfile(fid, dtype=dt,count=1)
    fid.seek(0,0)
    return f_dict[tagA[1]]

def acor_rd_m2( fid, frame_offset, ant_ind ):
    '''
    Output:
        Return the autocorrelation data of 16 channels of one antenna in Muser-I
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
        ant_ind: 0~43
    '''
    P = np.zeros(16);
    Arr = (np.arange(64)+1).reshape((64,1))
    Zarr = np.tile(np.array([1,1,1,1,0,0,0,0]).reshape(8,1),(8,1))
    Barr = np.repeat((np.arange(8)+1).reshape(1,8),8,1).T
    Sarr = np.tile(np.array([1,2,3,4]).reshape(4,1),(16,1))
    tablist = np.hstack([Arr, Zarr, Barr, Sarr]);
    dt = np.dtype('u4') #uint32
    frame_addr = frame_offset*204800; 
    auto_data_addr = 1888 + 12096; 
    zone_size = 12680; 
    block_size = 64; 
    slot_size = 8; 
    
    Nz = 1 if(tablist[ant_ind, 1]==1) else 2
    
    for kk in np.arange(Nz,16+1,2): # Nz:2:16:
        ant_auto_addr = frame_addr + auto_data_addr + (kk-1)*zone_size + \
                        (tablist[ant_ind, 2] - 1)*block_size + \
                         (tablist[ant_ind, 3] - 1)*slot_size
        ant_auto_addr = ant_auto_addr.astype(int)                        
        fid.seek(ant_auto_addr,0)                         
        if(Nz==1):
            P[kk-1] = np.fromfile(fid, count=1, dtype=dt)
            P[kk] = np.fromfile(fid, count=1, dtype=dt)
#            
#            P[kk-1] = fread(fid, 1, 'uint32')
#            P[kk] = fread(fid, 1, 'uint32')
        else:
            P[kk-2] = np.fromfile(fid, count=1, dtype=dt)
            P[kk-1] = np.fromfile(fid, count=1, dtype=dt)          
#            P[kk-2] = fread(fid, 1, 'uint32')
#            P[kk-1] = fread(fid, 1, 'uint32')          
    
    fid.seek(0,0)     
    return P                    
#    fseek(fid,0,'bof'); 

def xcor_rd_m2( fid, frame_offset, chan, antref, ant_list ):
    '''
    Output:
        Return the complex xcorrelation data of baseline: antref--ant_list[k]
    Input:
        fid: File handle of an opend Muser-I datafile.
        frame_offset: 0~19199, frame address in one datafile.
        chan: 0~15, channel number in 400 MHz
        antref: referenced antenna. 0~43
        antlist: [ant1, ant2, ant3...] 0~43
    '''
    zone_size = 12680;
    block_size = 12;
    slot_size = 2;
    dt = np.dtype('u1')
    xcor_I_out = np.zeros(len(ant_list))
    xcor_Q_out = np.zeros(len(ant_list))
    
#    ifelse=@(a,b,c)(a~=0)*b+(a==0)*c;
    xtable = np.zeros((64,64))
    xtable[0,:] = np.arange(64)
    for k in np.arange(63):
        xtable[k+1,:] = xtable[k,:] + 63-(k+1)
    
    for n in range(len(ant_list)):
        block_ind = np.ceil(xtable[antref,ant_list[n]]/2.)
        
        slot_ind = np.array([1,3,5]) if(np.mod(xtable[antref,ant_list[n]], 2)!=0) \
                                        else np.array([2,4,6])
#        slot_ind = ifelse(mod(xtable(antref,ant_list(n)), 2),[1 3 5],[2 4 6]);
    
        xcor_addr = frame_offset*204800 + 1888 + chan*zone_size \
                    + (block_ind - 1)*block_size + (slot_ind - 1)*slot_size    
        xcor_addr = xcor_addr.astype(int)                     
        fid.seek(0,0); fid.seek(xcor_addr[0],1)
        xcor_Q1 = np.fromfile(fid, count=1, dtype=dt) #[7:0]
        xcor_I1 = np.fromfile(fid, count=1, dtype=dt) #[7:0]
#        xcor_Q1 = fread(fid, 1, 'uint8') #[7:0]
#        xcor_I1 = fread(fid, 1, 'uint8')
        
        fid.seek(0,0);  fid.seek(xcor_addr[1],1)
        xcor_Q2 = np.fromfile(fid, count=1, dtype=dt) #[15:8]
        xcor_I2 = np.fromfile(fid, count=1, dtype=dt) 
#        xcor_Q2 = fread(fid, 1, 'uint8'); #[15:8]
#        xcor_I2 = fread(fid, 1, 'uint8');
        
        fid.seek(0,0);  fid.seek(xcor_addr[2],1)
        xcor_Q3 = np.fromfile(fid, count=1, dtype=dt) #[23:16]
        xcor_I3 = np.fromfile(fid, count=1, dtype=dt) 
#        xcor_Q3 = fread(fid, 1, 'uint8'); #[23:16]
#        xcor_I3 = fread(fid, 1, 'uint8');  
        
        xcor_I_out[n] = xcor_I3*256**2 + xcor_I2*256 + xcor_I1
        xcor_Q_out[n] = xcor_Q3*256**2 + xcor_Q2*256 + xcor_Q1
        xcor_I_out[n] = xcor_I_out[n]-2**24 if(xcor_I_out[n]>=2**23) else xcor_I_out[n]        
        xcor_Q_out[n] = xcor_Q_out[n]-2**24 if(xcor_Q_out[n]>=2**23) else xcor_Q_out[n]                
#        xcor_I_out[n] = ifelse(xcor_I_out(n)>2^23, xcor_I_out(n)-2^24, xcor_I_out(n));
#        xcor_Q_out(n) = ifelse(xcor_Q_out(n)>2^23, xcor_Q_out(n)-2^24, xcor_Q_out(n));
    
    fid.seek(0,0)
    return list(map(np.complex, xcor_I_out,xcor_Q_out))

def dly_rd_m1(fid, frame_offset):
    """
    Read out all MUSER-I antenna delay from a certain frame, including
    integer-part(ns) and fraction-part( 10^(-4) ns)    
    """
    N=40
    dly = np.zeros(N)
    BP = 100000*frame_offset + 99272 # address pointer, set to 99272 in each frame.
    bsize = 64
    dt = np.dtype('u1')
    for u in range(10):
        for v in range(4):
            BP_int = BP + 2*v + 8 + u*bsize # integer-part address
            BP_fra = BP + 2*v + u*bsize # fraction-part address
            fid.seek(BP_int,0)
            dly_int = np.fromfile(fid, count=2, dtype=dt)
            fid.seek(BP_fra,0)
            dly_fra = np.fromfile(fid, count=2, dtype=dt)            
            dly[4*u+v] = dly_int[1]*256 + dly_int[0] \
                        + (dly_fra[1]*256 + dly_fra[0])/10000.
    fid.seek(0,0)
    return dly

def dly_rd_m2(fid, frame_offset):    
    """
    Read out all MUSER-II antenna delay from a certain frame, including
    integer-part(ns) and fraction-part( 10^(-4) ns)    
    """
    return 0
    
#%%
#gps_str = '11010101001110101000101010101111'
#gps_st = {'yr':0,'mon':0}
#
#str1 = [gps_str[0:12],gps_str[12:16],gps_str[16:21],gps_str[21:26],gps_str[26:32],\
#        gps_str[32:38],gps_str[38:48],gps_str[48:58],gps_str[58:64]]
#
#gps_st['yr'], gps_st['mon'],gps_st['day'],gps_st['hr'],gps_st['min'], \
#gps_st['sec'],gps_st['msec'],gps_st['usec'],gps_st['nsec']  \
#= map(int,str1,[2,2,2,2,2,2,2,2,2])
#gps_st['yr']=gps_st['yr']+2000

#t1 = {'yr':2016,'mon':8,'day':2,'hr':3,'min':7,'sec':9,'msec':6,'usec':12,'nsec':9}
#gps_isot = '{yr}-{mon:02}-{day:02}T{hr:02}:{min:02}:{sec:02}.{msec:03}{usec:03}'.format(**t1)
#print gps_isot1
##%%
#gps_isot = '{:0{4}}-{:0{2}}-{:0{2}}T{:0{2}}:{:0{2}}:{:0{2}}.{:0{3}}{:0{3}}{:0{3}}'\
#            .format(**gps_st)    