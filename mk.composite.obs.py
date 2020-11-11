#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import util
from bisect import bisect_left
from detect_fsub import *
import d4PDF
import JRA55
import IBTrACS
from myfunc.regrid import Regrid
import myfunc.grids as grids
#--------------------------------------

iYM = [1995,5]
eYM = [2018,12]
#iYM = [1980,1]
#eYM = [2018,12]

lYM = util.ret_lYM(iYM,eYM)
#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------

figdir  = '/home/utsumi/temp/bams2020'
compbasedir= '/home/utsumi/temp/bams2020/composite'
jrabaseDir = '/home/utsumi/mnt/lab_data2/JRA55/'
#----------------------------------

a1latin = JRA55.Lat125(crd='sa')
a1lonin = JRA55.Lon125(crd='sa')

a1latout = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonout = d4PDF.Lon()   # dlon = 0.5625  # 0, 0.5625, ..

latout0 = a1latout[0]
lonout0 = a1lonout[0]
nyout   = len(a1latout)
nxout   = len(a1lonout)
dlatout = 0.5615674
dlonout = 0.5625

jra_anl_p    = JRA55.anl_p125(dbbaseDir=jrabaseDir)
jra_fcst_phy2m = JRA55.fcst_phy2m125(dbbaseDir=jrabaseDir)
jra_anl_surf = JRA55.anl_surf125(dbbaseDir=jrabaseDir)

nyin    = len(a1latin)
nxin    = len(a1lonin)

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

    a3slp  = []
    a3prec = []
    a3wind = []
    a3wind500=[]
    a3wind850=[]
    a3t500 = []
    a3t850 = []
    a1vort = []
    a1tsum = []
    a1lat  = []
    a1lon  = []
    a1time = [] 

    lDTimeBst = np.sort(dlonlat.keys())
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

            x = int((lon-lonout0)/dlonout)
            y = int((lat-latout0)/dlatout)
            lx.append(x)
            ly.append(y)
            llat.append(lat)
            llon.append(lon)

        a1x = np.array(lx)
        a1y = np.array(ly)

        if len(a1x)==0: continue

        #--- Loading data ----------
        print DTime

        a2vort = jra_anl_p.load_6hr(vname='relv', DTime=DTime, plev=850, crd='sa')
        a2t500  = jra_anl_p.load_6hr(vname='tmp', DTime=DTime, plev=500, crd='sa') 
        a2t850  = jra_anl_p.load_6hr(vname='tmp', DTime=DTime, plev=850, crd='sa') 
        a2u850  = jra_anl_p.load_6hr(vname='ugrd', DTime=DTime, plev=850, crd='sa') 
        a2v850  = jra_anl_p.load_6hr(vname='vgrd', DTime=DTime, plev=850, crd='sa') 
        a2u500  = jra_anl_p.load_6hr(vname='ugrd', DTime=DTime, plev=500, crd='sa') 
        a2v500  = jra_anl_p.load_6hr(vname='vgrd', DTime=DTime, plev=500, crd='sa') 

        a2u    = jra_anl_surf.load_6hr(vname='10 metre U wind component', DTime=DTime, crd='sa')
        a2v    = jra_anl_surf.load_6hr(vname='10 metre V wind component', DTime=DTime, crd='sa')
        a2slp  = jra_anl_surf.load_6hr(vname='Mean sea level pressure', DTime=DTime, crd='sa')

        a2prec0 = jra_fcst_phy2m.load_3hr(vname='Mean total precipitation', DTime=DTime, fcst='forward', crd='sa')
        a2prec1 = jra_fcst_phy2m.load_3hr(vname='Mean total precipitation', DTime=DTime, fcst='backward', crd='sa')
        a2prec  = 0.5*(a2prec0 + a2prec1) / 24.  # mm/day --> mm/hour
        #-- Downscaling -----------
        a2slp  = Regrid.biIntp(a1latin, a1lonin, a2slp, a1latout, a1lonout)[0]
        a2vort = Regrid.biIntp(a1latin, a1lonin, a2vort, a1latout, a1lonout)[0]
        a2u    = Regrid.biIntp(a1latin, a1lonin, a2u, a1latout, a1lonout)[0]
        a2v    = Regrid.biIntp(a1latin, a1lonin, a2v, a1latout, a1lonout)[0]
        a2prec = Regrid.biIntp(a1latin, a1lonin, a2prec, a1latout, a1lonout)[0]
        a2u500 = Regrid.biIntp(a1latin, a1lonin, a2u500, a1latout, a1lonout)[0]
        a2v500 = Regrid.biIntp(a1latin, a1lonin, a2v500, a1latout, a1lonout)[0]
        a2u850 = Regrid.biIntp(a1latin, a1lonin, a2u850, a1latout, a1lonout)[0]
        a2v850 = Regrid.biIntp(a1latin, a1lonin, a2v850, a1latout, a1lonout)[0]
        a2t500 = Regrid.biIntp(a1latin, a1lonin, a2t500, a1latout, a1lonout)[0]
        a2t850 = Regrid.biIntp(a1latin, a1lonin, a2t850, a1latout, a1lonout)[0]
        #-- Wind speed ------------
        a2wind = np.sqrt(np.square(a2u) + np.square(a2v))
        a2wind500=np.sqrt(np.square(a2u500) + np.square(a2v500))
        a2wind850=np.sqrt(np.square(a2u850) + np.square(a2v850))

        #-- Max search rvort ------
        a2vort[:nyout/2] = -a2vort[:nyout/2] 
        a2vort = grids.karnel_pooling_map2D_global(a2vort, dy=1, dx=1, func='max')

        #-- 1D array -------------
        a1vort = a1vort + a2vort[a1y,a1x].tolist()
        a1lat  = a1lat  + llat
        a1lon  = a1lon  + llon
 
        for (y,x) in zip(a1y,a1x): 

            iy = y-dgridy
            ey = y+dgridy
            ix = x-dgridx
            ex = x+dgridx
    
            a2slptmp  = a2slp[iy:ey+1,ix:ex+1]
            a2prectmp = a2prec[iy:ey+1,ix:ex+1]
            a2windtmp = a2wind[iy:ey+1,ix:ex+1]
            a2wind500tmp = a2wind500[iy:ey+1,ix:ex+1]
            a2wind850tmp = a2wind850[iy:ey+1,ix:ex+1]
            a2t850tmp = a2t850[iy:ey+1,ix:ex+1]
            a2t500tmp = a2t500[iy:ey+1,ix:ex+1]

            a3slp.append(a2slptmp)
            a3prec.append(a2prectmp)
            a3wind.append(a2windtmp)
            a3wind500.append(a2wind500tmp)
            a3wind850.append(a2wind850tmp)
            a3t500.append(a2t500tmp)
            a3t850.append(a2t850tmp)
            a1time.append(DTime) 

    dout = {}
    dout['slp']= np.array(a3slp).astype('float32')
    dout['prec']= np.array(a3prec).astype('float32')
    dout['wind']= np.array(a3wind).astype('float32')
    dout['wind500']= np.array(a3wind500).astype('float32')
    dout['wind850']= np.array(a3wind850).astype('float32')
    dout['t500']= np.array(a3t500).astype('float32')
    dout['t850']= np.array(a3t850).astype('float32')
    dout['vort']= np.array(a1vort).astype('float32')
    dout['lat' ]= np.array(a1lat).astype('float32')
    dout['lon' ]= np.array(a1lon).astype('float32')
    dout['time']= np.array(a1time)

    #--- Save file ---------
    compdir = compbasedir + '/obs-tc'
    util.mk_dir(compdir)

    comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)
    with open(comppath,'wb') as f:
         
        pickle.dump(dout, f)
        print comppath
         
