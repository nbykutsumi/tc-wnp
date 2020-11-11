# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
from myfunc.regrid import Regrid
import myfunc.grids as grids
#--------------------------------------
#iYM = [1980,1]
iYM = [1986,10]
eYM = [2018,12]

#iYM = [1980,1]
#eYM = [1980,1]

lYM = util.ret_lYM(iYM,eYM)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "JRA55"
model   = "__"
run     = "__"
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}

detectName = 'wsd_d4pdf_20200428'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))
JRA55       = import_module("%s.JRA55"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))

wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/JRA55'
figdir  = '/home/utsumi/temp/bams2020'
compbasedir= '/home/utsumi/temp/bams2020/composite'
jrabaseDir = '/home/utsumi/mnt/lab_data2/JRA55/'

#----------------------------------
a1latin = JRA55.Lat125(crd='sa')   # dlat = 1.25
a1lonin = JRA55.Lon125(crd='sa')   # dlon = 1.25
a1latout = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonout = d4PDF.Lon()   # dlon = 0.5625  # 0, 0.5625, ..

nyin    = len(a1latin)
nxin    = len(a1lonin)
nyout   = len(a1latout)
nxout   = len(a1lonout)

miss =-9999
miss_int= -9999

jra_anl_p    = JRA55.anl_p125(dbbaseDir=jrabaseDir)
jra_fcst_phy2m = JRA55.fcst_phy2m125(dbbaseDir=jrabaseDir)
jra_anl_surf = JRA55.anl_surf125(dbbaseDir=jrabaseDir)


#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

#************************************
# JRA55 objective detection for ExC
#************************************
const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = a1latin
const['Lon'] = a1lonin

thsst   = const['thsst']  - 273.15   # degC
exrvort  = const['exrvort'] 
thdura  = const['thdura']  
exrvortout = const['exrvort']*1.0e+5

slabel = 'sst-%d.ex-%.2f.1f-du-%02d'%(thsst, exrvortout, thdura)

wsDir= wsbaseDir + '/%s'%(run)

cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)

for YM in lYM:
    Year,Mon = YM
    lDTime = util.ret_lDTime_fromYM([Year,Mon],[Year,Mon], timedelta(hours=6), hour0=0)

    #--- load TC dictionary ---
    dcxy, _ = cy.mkInstDictC_bstTC(YM,YM,varname='vortlw')

    #--- load TC ipos, idate, dura ------
    dcxy_ipos , _ = cy.mkInstDictC_bstTC(YM,YM,varname='ipos')
    dcxy_idate, _ = cy.mkInstDictC_bstTC(YM,YM,varname='idate')
    dcxy_dura , _ = cy.mkInstDictC_bstTC(YM,YM,varname='dura')
    #--------------------------    

    a3prec = []
    a3wind = []
    a3wind500=[]
    a3wind850=[]
    a3t500 = []
    a3t850 = []
    a3slp  = []
    a1lat  = []
    a1lon  = []
    a1time = [] 
    a1vort = [] 
    a1ipos = []
    a1idate= []
    a1dura = []

    #lDTime = [datetime(1980,1,5,0)]  # test
    #lDTime = util.ret_lDTime(datetime(1980,1,1,0), datetime(1980,1,5,0),timedelta(hours=6))  # test
    for DTime in lDTime:
        #--------------------------
        lcxy      = dcxy[DTime]
        lcxy_ipos = dcxy_ipos[DTime] 
        lcxy_idate= dcxy_idate[DTime] 
        lcxy_dura = dcxy_dura [DTime] 

        if len(lcxy)==0: continue

        print 'org=',len(lcxy)
        lcx, lcy, lcvar = map(np.array,zip(*lcxy))
        _,_, lipos = map(np.array, zip(*lcxy_ipos))
        _,_, lidate= map(np.array, zip(*lcxy_idate))
        _,_, ldura = map(np.array, zip(*lcxy_dura ))


        llattmp = a1latin[lcy] 
        llontmp = a1lonin[lcx]

        a1flaglat= ma.masked_inside(llattmp,lllat,urlat).mask
        a1flaglon= ma.masked_inside(llontmp,lllon,urlon).mask
        a1flag   = a1flaglat * a1flaglon

        if a1flag.sum()==0: continue

        llattmp = llattmp[a1flag]
        llontmp = llontmp[a1flag]
        lytmp   = lcy[a1flag]
        lxtmp   = lcx[a1flag]
        lvorttmp= lcvar[a1flag]
        lipostmp= lipos[a1flag]
        lidatetmp= lidate[a1flag]
        lduratmp = ldura[a1flag]

        print len(llattmp)
        a1lat  = a1lat + llattmp.tolist() 
        a1lon  = a1lon + llontmp.tolist()
        a1vort = a1vort+ lvorttmp.tolist()
        a1ipos = a1ipos+ lipostmp.tolist()
        a1idate= a1idate+ lidatetmp.tolist()
        a1dura = a1dura + lduratmp.tolist()

        #--- Loading JRA55 data ----------
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

        for i in range(len(lxtmp)):
            lat = llattmp[i]
            lon = llontmp[i]

            #yout = bisect_left(a1latout, lat)
            yout = int((lat - a1latout[0] + 0.5615674*0.5)/0.5615674)

            xout = int((lon - a1lonout[0] + 0.5625*0.5)/0.5625)

            iy = yout-dgridy
            ey = yout+dgridy
            ix = xout-dgridx
            ex = xout+dgridx

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
    compdir = compbasedir + '/obs-exc'
    util.mk_dir(compdir)

    comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)
    with open(comppath,'wb') as f:

        pickle.dump(dout, f)
        print comppath



# %%
