# %%
import numpy as np
import d4PDF, JRA55
from numpy import *
from datetime import datetime, timedelta
import myfunc.util  as util
import myfunc.grids as grids
import os, sys
from collections import deque
import pickle, calendar
import bisect

#obscalc = True
obscalc = False
#simcalc = True
simcalc = False

iYM = [2000,1]
eYM = [2002,12]
#eYM = [2010,1]
lYM = util.ret_lYM(iYM,eYM)
BBox = [[5,100],[60,180]]
[[lllat,lllon],[urlat,urlon]] = BBox
picklebasedir = '/home/utsumi/temp/bams2020/pickle'
csvdir        = '/home/utsumi/temp/bams2020'
vname = 'rvort'

#*******************************************
# Function
#-------------------------------------------
def calc_distance():
    a1lat = d4PDF.Lat()
    a1lon = d4PDF.Lon()

    ny = len(a1lat)
    nx = len(a1lon)

    RADEARTH = 6371
    DTR = 0.017453

    lon1, lon2 = a1lon[0], a1lon[2]
    lat1, lat2 = a1lat, a1lat
    a1distx = RADEARTH*np.arccos(np.cos(DTR*lon1-DTR*lon2)*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))

    lon1, lon2 = 0, 0
    lat1, lat2 = a1lat[79], a1lat[81]
    disty = RADEARTH*np.arccos(np.cos(DTR*lon1-DTR*lon2)*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))
    a1disty = np.full(ny, disty)

    a2disty = np.repeat(a1disty.reshape(-1,1), nx, axis=1)
    a2distx = np.repeat(a1distx.reshape(-1,1), nx, axis=1)

    return a2disty, a2distx

a2disty, a2distx = calc_distance()


def calc_a2rvort(DTime, scen, ens):
    a2u = d4atm.load_6hr(vname='U850', scen=scen, ens=ens, DTime=DTime)  # missing masked
    a2v = d4atm.load_6hr(vname='V850', scen=scen, ens=ens, DTime=DTime)  # missing masked

    a2un = ma.masked_equal( grids.shift_map_global(a2u, dy=-1, dx=0, miss=-9999.) ,-9999.)
    a2us = ma.masked_equal( grids.shift_map_global(a2u, dy=1, dx=0, miss=-9999.)  ,-9999.)
    a2ve = ma.masked_equal( grids.shift_map_global(a2v, dy=0, dx=-1, miss=-9999.) ,-9999.)
    a2vw = ma.masked_equal( grids.shift_map_global(a2v, dy=0, dx=1, miss=-9999.)  ,-9999.)

    a2rvort = ((a2ve-a2vw)/a2distx - (a2un-a2us)/a2disty)*1e-3

    return a2rvort




#*************************************
# JRA55
#-------------------------------------
vnamejra = {'rvort':'relv'}[vname]
jrabaseDir = '/home/utsumi/mnt/lab_data2/JRA55/'
jra = JRA55.anl_p125(dbbaseDir=jrabaseDir)
a1latjra = JRA55.Lat125(crd='sa')
a1lonjra = JRA55.Lon125(crd='sa')
a2lonjra, a2latjra = np.meshgrid(a1lonjra, a1latjra)
dlat,dlon= JRA55.dlatlon125()
latjra0  = a1latjra[0]
lonjra0  = a1lonjra[0]
nyjra    = len(a1latjra)
nxjra    = len(a1lonjra)
iy = int((lllat - latjra0)/dlat)
ix = int((lllon - lonjra0)/dlon)
ey = int((urlat - latjra0)/dlat)
ex = int((urlon - lonjra0)/dlon)

avar = deque([])
alat = deque([])
alon = deque([])
for (Year,Mon) in lYM:
    if obscalc is not True: continue

    eday = calendar.monthrange(Year,Mon)[1] 
    iDTimejra = datetime(Year,Mon,1)
    eDTimejra = datetime(Year,Mon,eday)
    lDTimejra = util.ret_lDTime(iDTimejra,eDTimejra,timedelta(hours=6))

    for DTime in lDTimejra:
        print 'JRA55',DTime
        a2var = jra.load_6hr(vname=vnamejra, DTime=DTime, plev=850, crd='sa')
        a2var[:nyjra/2] = (-ma.masked_equal(a2var[:nyjra/2], -9999.)).filled(-9999.) # Flip signs of southern hemisphere
        a2vartmp = a2var[iy:ey+1, ix:ex+1]
        a2lattmp = a2latjra[iy:ey+1, ix:ex+1]
        a2lontmp = a2lonjra[iy:ey+1, ix:ex+1]
        
        a2flag   = ma.masked_greater(a2vartmp,0).mask
        avar.extend(a2vartmp[a2flag])
        alat.extend(a2lattmp[a2flag])
        alon.extend(a2lontmp[a2flag])


obsdir = picklebasedir + '/pdfmatch-%s'%(vname)
util.mk_dir(obsdir)
obspath= obsdir + '/jra55.pickle'

if obscalc is True:
    dout = {}
    dout['BBox'] = BBox
    dout['lYM']  = lYM
    dout['a1var'] = np.array(avar).astype('float32')
    dout['a1lat'] = np.array(alat).astype('float32')
    dout['a1lon'] = np.array(alon).astype('float32')
    with open(obspath, 'wb') as f:
        pickle.dump(dout, f)
    print obspath

dout = {}

#************************************************
# d4PDF
#------------------------------------------------
prj    = 'd4PDF'
model  = '__'
expr   = 'XX'
scen   = 'HPB'
ens    = 1
res    = '320x640'

wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
d4atm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir='/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM')

a1latsim = d4PDF.Lat()
a1lonsim = d4PDF.Lon()
a2lonsim, a2latsim = np.meshgrid(a1lonsim, a1latsim)
nysim    = len(a1latsim)
nxsim    = len(a1lonsim)
iy = bisect.bisect_left(a1latsim, lllat)
ey = bisect.bisect_left(a1latsim, urlat)
ix = bisect.bisect_left(a1lonsim, lllon)
ex = bisect.bisect_left(a1lonsim, urlon)

avar = deque([])
alat = deque([])
alon = deque([])
for (Year,Mon) in lYM:
    if simcalc is not True: continue

    eday = calendar.monthrange(Year,Mon)[1] 
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eday)
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=6))

    for DTime in lDTime:
        print 'd4PDF', DTime
        if vname == 'rvort':
            a2var    = calc_a2rvort(DTime, scen, ens)
            a2var[:160] = -ma.masked_equal(a2var[:160],-9999.)
        else:
            print 'check vname',vname
            sys.exit()
        a2vartmp = a2var[iy:ey+1, ix:ex+1]
        a2lattmp = a2latsim[iy:ey+1, ix:ex+1]
        a2lontmp = a2lonsim[iy:ey+1, ix:ex+1]
        
        a2flag   = ma.masked_greater(a2vartmp,0).mask
        avar.extend(a2vartmp[a2flag])
        alat.extend(a2lattmp[a2flag])
        alon.extend(a2lontmp[a2flag])

simdir = picklebasedir + '/pdfmatch-%s'%(vname)
util.mk_dir(simdir)
simpath= simdir + '/d4PDF.pickle'

if simcalc is True:
    dout = {}
    dout['BBox'] = BBox
    dout['lYM']  = lYM
    dout['a1var'] = np.array(avar).astype('float32')
    dout['a1lat'] = np.array(alat).astype('float32')
    dout['a1lon'] = np.array(alon).astype('float32')
    with open(simpath, 'wb') as f:
        pickle.dump(dout, f)
    print simpath


#*********************************************************
# Figure
#---------------------------------------------------------
#lpercent = np.arange(99.4,100,0.0001)
lpercent = list(np.arange(50,90,0.1)) + list(np.arange(90,99,0.01)) + list(np.arange(99,99.9,0.001)) + list(np.arange(99.9,100,0.0001))
with open(obspath, 'rb') as f:
    ddat = pickle.load(f)
    avar = ddat['a1var']
aobs = np.percentile(avar, lpercent)
print 'len obs=',len(avar)

with open(simpath, 'rb') as f:
    ddat = pickle.load(f)
    avar = ma.masked_invalid(ddat['a1var']).compressed()
asim = np.percentile(avar, lpercent)
print 'len sim=',len(avar)

ddat = None
avar = None
lout = [['percentile,obs,sim']]
for i,percent in enumerate(lpercent):
    lout.append([percent, aobs[i],asim[i]])

sout = util.list2csv(lout)
csvpath = csvdir + '/%s.percentile.csv'%(vname)
f=open(csvpath,'w'); f.write(sout); f.close()
print csvpath
# %%
