#----------------------------------
import sys, os, socket
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import util
from bisect import bisect_left
from detect_fsub import *
from netCDF4 import Dataset
import IBTrACS
#--------------------------------------

iYM = [1980,1]
eYM = [2002,12]
#eYM = [2018,12]
#iYM = [1980,1]
#eYM = [2018,12]

lYM = util.ret_lYM(iYM,eYM)
#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------

figdir  = '/home/utsumi/temp/bams2020'
compbasedir= '/home/utsumi/temp/bams2020/composite'

myhost=socket.gethostname()
if myhost =='shui':
    erabaseDir =  '/tank/utsumi/era5/deg05625'
elif myhost=='well':
    erabaseDir = '/media/disk2/data/era5/deg05625'
else:
    print('check myhost',myhost)
    sys.exit()

def read_ave_2d_hour_trange(vname,iDTime,eDTime, mode="3d"):
    """ 
        total and average: BACKWARD ending at the time stamp
        This function does not intend efficient calculation.
    """
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=1))

    a3out = []
    for DTime in lDTime:
        Year,Mon,Day,Hour = DTime.timetuple()[:4]
        srcDir = erabaseDir + '/%s/%04d%02d'%(vname,Year,Mon)
        srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(vname,Year,Mon,Day)
        with Dataset(srcPath,'r') as np:
            a3out.append(np.variables[vname][Hour])

    if mode=="3d":
        return ma.asarray(a3out)
    elif mode=="ave":
        return ma.asarray(a3out).mean(axis=0)
    elif mode=="sum":
        return ma.asarray(a3out).sum(axis=0)
    else:
        print("check mode",mode)
        sys.exit()


def read_var_2d_hour(vname,DTime):
    """ 
        total and average: BACKWARD ending at the time stamp
        This function does not intend efficient calculation.
    """
    Year,Mon,Day,Hour = DTime.timetuple()[:4]
    srcDir = erabaseDir + '/%s/%04d%02d'%(vname,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(vname,Year,Mon,Day)
    with Dataset(srcPath,'r') as np:
        a2var = np.variables[vname][Hour]
    return a2var

def read_fix(vname):
    srcDir = erabaseDir + '/fix'
    srcPath= srcDir + '/%s.nc'%(vname,Year,Mon,Day)
    with Dataset(srcPath,'r') as np:
        a2var = np.variables[vname][0]
    return a2var


#----------------------------------

latin0  = -90  # after flipping
lonin0  = 0
nyin    = 321
nxin    = 640
dlat    = 0.5625
dlon    = 0.5625


miss =-9999
miss_int= -9999

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]
#[[lllat,lllon],[urlat,urlon]] = [[20,120],[47,150]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]
#**************************
#Read best track
#---------------------------
for YM in lYM:
    Year,Mon = YM
    lDTime = util.ret_lDTime_fromYM([Year,Mon],[Year,Mon], timedelta(hours=6), hour0=0)
    iDTime = lDTime[0]
    eDTime = lDTime[-1]

    dbst   = {}
    bst    = IBTrACS.IBTrACS()
    dlonlat = bst.ret_dlonlat(iDTime,eDTime)
    #dlonlat = bst.ret_dlonlat(iDTime,eDTime,natures=['TC','NR'])

    #a3slp  = []
    a3prec0= []
    a3prec1= []
    a3prec2= []
    a3prec3= []
    a3wind = []
    a3land = []
    #a1vort = []
    #a1tsum = []
    a1lat  = []
    a1lon  = []
    a1time = [] 

    lDTimeBst = np.sort(list(dlonlat.keys()))
    for DTime in lDTimeBst:
        ltclonlat = dlonlat[DTime]
        #--- Lat, Lon screening ----
        lx = []
        ly = []
        llat=[]
        llon=[]
        for (lon,lat) in ltclonlat: 
            if (lon<lllon)or(urlon<lon):
                continue
            if (lat<lllat)or(urlat<lat):
                continue   

            x = int((lon-lonin0)/dlon)
            y = int((lat-latin0)/dlat)
            lx.append(x)
            ly.append(y)
            llat.append(lat)
            llon.append(lon)

        a1x = np.array(lx)
        a1y = np.array(ly)

        if len(a1x)==0: continue

        #--- Loading data ----------
        print(DTime)

        #a2vort = jra_anl_p.load_6hr(vname='relv', DTime=DTime, plev=850, crd='sa')
        a2u    = np.flipud(read_var_2d_hour(vname='u10', DTime=DTime))
        a2v    = np.flipud(read_var_2d_hour(vname='v10', DTime=DTime))


        a2prec0= np.flipud(read_ave_2d_hour_trange(vname='tp', iDTime=DTime+timedelta(hours=-11), eDTime=DTime+timedelta(hours=-6),mode="ave"))  # backward accumulation ending at the time stamp, m/h
        a2prec1= np.flipud(read_ave_2d_hour_trange(vname='tp', iDTime=DTime+timedelta(hours=-5),  eDTime=DTime+timedelta(hours=0),mode="ave"))  # backward accumulation ending at the time stamp, m/h
        a2prec2= np.flipud(read_ave_2d_hour_trange(vname='tp', iDTime=DTime+timedelta(hours=+1),  eDTime=DTime+timedelta(hours=6),mode="ave"))  # backward accumulation ending at the time stamp, m/h
        a2prec3= np.flipud(read_ave_2d_hour_trange(vname='tp', iDTime=DTime+timedelta(hours=+7),  eDTime=DTime+timedelta(hours=12),mode="ave"))  # backward accumulation ending at the time stamp, m/h


        #-- Wind speed ------------
        a2wind = np.sqrt(np.square(a2u) + np.square(a2v))

        ##-- Max search rvort ------
        #a2vort[:nyout/2] = -a2vort[:nyout/2] 
        #a2vort = grids.karnel_pooling_map2D_global(a2vort, dy=1, dx=1, func='max')

        #-- 1D array -------------
        #a1vort = a1vort + a2vort[a1y,a1x].tolist()
        a1lat  = a1lat  + llat
        a1lon  = a1lon  + llon
 
        for (y,x) in zip(a1y,a1x): 

            iy = y-dgridy
            ey = y+dgridy
            ix = x-dgridx
            ex = x+dgridx
    
            #a2slptmp  = a2slp[iy:ey+1,ix:ex+1]
            a2prectmp0 = a2prec0[iy:ey+1,ix:ex+1]
            a2prectmp1 = a2prec1[iy:ey+1,ix:ex+1]
            a2prectmp2 = a2prec2[iy:ey+1,ix:ex+1]
            a2prectmp3 = a2prec3[iy:ey+1,ix:ex+1]
            a2windtmp = a2wind[iy:ey+1,ix:ex+1]

            #a3slp.append(a2slptmp)
            a3prec0.append(a2prectmp0)
            a3prec1.append(a2prectmp1)
            a3prec2.append(a2prectmp2)
            a3prec3.append(a2prectmp3)

            a3wind.append(a2windtmp)
            a1time.append(DTime) 

    dout = {}
    #dout['slp']= np.array(a3slp).astype('float32')
    dout['prec-12']= np.array(a3prec0).astype('float32')   # meter/h
    dout['prec-06']= np.array(a3prec1).astype('float32')   # meter/h
    dout['prec000']= np.array(a3prec2).astype('float32')   # meter/h
    dout['prec006']= np.array(a3prec3).astype('float32')   # meter/h
    dout['wind']= np.array(a3wind).astype('float32')
    dout['land']= np.array(a3land).astype('float32')
    #dout['vort']= np.array(a1vort).astype('float32')
    dout['lat' ]= np.array(a1lat).astype('float32')
    dout['lon' ]= np.array(a1lon).astype('float32')
    dout['time']= np.array(a1time)

    #--- Save file ---------
    compdir = compbasedir + '/obs-era5-tc/%04d'%(Year)
    util.mk_dir(compdir)
    for vname in list(dout.keys()):
        opath = compdir + '/%s.%04d.%02d.npy'%(vname,Year,Mon)
        np.save(opath, dout[vname])
        print(opath)
         
