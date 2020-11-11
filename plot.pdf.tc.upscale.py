# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import myfunc.util as util
import myfunc.grids as grids
import Cyclone
import ConstCyclone
import d4PDF
import IBTrACS, JRA55
from datetime import datetime, timedelta
import sys, os
import pickle
#************************************************
#obsflag = True
obsflag = False
#simflag = True
simflag = False

d4atm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir='/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM')

#************************************************
# Functions
#------------------------------------------------
def calc_upscale_distance():
    a1latorg = d4PDF.Lat()
    a1lonorg = d4PDF.Lon()
    a1latup  = a1latorg.reshape(-1,2).mean(axis=1)
    a1lonup  = a1lonorg.reshape(-1,2).mean(axis=1)

    nyup = len(a1latup)
    nxup = len(a1lonup)

    RADEARTH = 6371
    DTR = 0.017453

    lon1, lon2 = a1lonup[0], a1lonup[2]
    lat1, lat2 = a1latup, a1latup
    a1distx = RADEARTH*np.arccos(np.cos(DTR*lon1-DTR*lon2)*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))

    lon1, lon2 = 0, 0
    lat1, lat2 = a1latup[79], a1latup[81]
    disty = RADEARTH*np.arccos(np.cos(DTR*lon1-DTR*lon2)*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))
    a1disty = np.full(nyup, disty)

    a2disty = np.repeat(a1disty.reshape(-1,1), nxup, axis=1)
    a2distx = np.repeat(a1distx.reshape(-1,1), nxup, axis=1)

    return a2disty, a2distx

a2distupy, a2distupx = calc_upscale_distance()

def calc_upscale_a2rvort(DTime, scen, ens):
    a2uorg = d4atm.load_6hr(vname='U850', scen=scen, ens=ens, DTime=DTime)  # missing masked
    a2vorg = d4atm.load_6hr(vname='V850', scen=scen, ens=ens, DTime=DTime)  # missing masked

    a2uup = grids.fold2d(a2uorg, [160,320]).mean(axis=(2,3))  # keep masks
    a2vup = grids.fold2d(a2vorg, [160,320]).mean(axis=(2,3))  # keep masks

    a2un = ma.masked_equal( grids.shift_map_global(a2uup, dy=-1, dx=0, miss=-9999.) ,-9999.)
    a2us = ma.masked_equal( grids.shift_map_global(a2uup, dy=1, dx=0, miss=-9999.)  ,-9999.)
    a2ve = ma.masked_equal( grids.shift_map_global(a2vup, dy=0, dx=-1, miss=-9999.) ,-9999.)
    a2vw = ma.masked_equal( grids.shift_map_global(a2vup, dy=0, dx=1, miss=-9999.)  ,-9999.)

    a2rvortup = ((a2ve-a2vw)/a2distupx - (a2un-a2us)/a2distupy)*1e-3


    return a2rvortup
#************************************************
# IBTrACS
#------------------------------------------------
iDTime = datetime(1980,1,1,0)
#eDTime = datetime(1981,1,1,0)
eDTime = datetime(2019,12,31,18)
bst = IBTrACS.IBTrACS()
bstlonlat = bst.ret_dlonlat(iDTime, eDTime)
lDTimeBst = np.sort(bstlonlat.keys())


jrabaseDir = '/home/utsumi/mnt/lab_data2/JRA55/'
jra = JRA55.anl_p125(dbbaseDir=jrabaseDir)
vname = 'relv'
a1latjra = JRA55.Lat125(crd='sa')
a1lonjra = JRA55.Lon125(crd='sa')
dlat,dlon= JRA55.dlatlon125()
latjra0  = a1latjra[0]
lonjra0  = a1lonjra[0]

a1varout = []
a1latout = []
a1lonout = []

if obsflag == True:
    for DTime in lDTimeBst:
        print DTime
        nyjra = len(JRA55.Lat125())
        avar = jra.load_6hr(vname=vname, DTime=DTime, plev=850, crd='sa')
        avar[:nyjra/2] = (-ma.masked_equal(avar[:nyjra/2], -9999.)).filled(-9999.) # Flip signs of southern hemisphere
        avar = grids.karnel_pooling_map2D_global(avar, dy=1,dx=1,func='max',miss_in=-9999.)  # Find maximum in 3x3 grid boxes
        llonlat = bstlonlat[DTime]
        for (lon,lat) in llonlat:
            if lon<0: lon=lon+360
            y = int((lat - latjra0)/dlat)
            x = int((lon - lonjra0)/dlon)

            a1varout.append(avar[y,x])

            a1latout.append(lat)
            a1lonout.append(lon)

    dobs = {}
    dobs[vname] = np.array(a1varout)
    dobs['Lat'] = np.array(a1latout)
    dobs['Lon'] = np.array(a1lonout)

#plt.hist(a1varout)
#plt.show()
#sys.exit()
#-- Save -------
pickledir = '/home/utsumi/temp/bams2020/pickle'
util.mk_dir(pickledir)
obspath = pickledir + '/obs.%s.pickle'%(vname)
if obsflag==True:
    with open(obspath, 'wb') as f:
        pickle.dump(dobs, f)
    print obspath

#************************************************
# d4PDF
#------------------------------------------------
prj    = 'd4PDF'
model  = '__'
expr   = 'XX'
scen   = 'HPB'
lens   = range(1,9+1)
#lens   = [1]
res    = '320x640'

iYM = [2000,8]
eYM = [2010,12]
#eYM = [2000,9]
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

const = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = d4PDF.Lat()
const['Lon'] = d4PDF.Lon()

thsst   = 27  # degC
exrvort = 3 * 1e-5
tcrvort = 4 * 1e-5
thwcore = 0
thdura  = 60
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5
const['thsst']   = thsst + 273.15
const['exrvort'] = exrvort
const['tcrvort'] = tcrvort
const['thwcore'] = thwcore
const['thdura']  = thdura
#--------------------------------------------
slabel = 'upscale.sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thdura)

simdir = '/home/utsumi/temp/bams2020/pickle/%s'%(slabel)
util.mk_dir(simdir)

vname = 'vortlw'
for ens in lens:
    if simflag is not True:
        continue

    print 'calc d4PDF',ens
    wsDir = wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    cy = Cyclone.Cyclone(baseDir=wsDir, const=const)
    dexcxy, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname=vname)

    #-- Find rvort from upscaled d4PDF --
    ldtimetc = dtcxy.keys()
    llxyv    = dtcxy.values()
    a1var    = []
    a1x      = []
    a1y      = []
    for itc,dtimetc in enumerate(ldtimetc):
        lxyv = llxyv[itc]
        if len(lxyv)==0: continue

        lxy = np.array(llxyv[itc])[:,:2].astype('int')

        a2rvortup= calc_upscale_a2rvort(dtimetc, scen, ens)
        a2rvortup[:80] = -ma.masked_equal(a2rvortup[:80],-9999.)

        a2rvortup = grids.karnel_pooling_map2D_global(a2rvortup, dy=1, dx=1, func='max', miss_in=-9999., miss_out=-9999.)

        for (x,y) in lxy:
            yup = int(y/2)
            xup = int(x/2)
            a1var.append(a2rvortup[yup,xup])
            a1x.append(x)
            a1y.append(y)

            
    #---------------
    a1latsim = d4PDF.Lat()[a1y]
    a1lonsim = d4PDF.Lon()[a1x]

    dsim = {}
    dsim[vname] = np.array(a1var)
    dsim['Lat'] = np.array(a1latsim)
    dsim['Lon'] = np.array(a1lonsim)
    dsim['iYM'] = iYM
    dsim['eYM'] = eYM

    simpath = simdir + '/%s.%s.%03d.pickle'%(vname,scen,ens)

    with open(simpath, 'wb') as f:
        pickle.dump(dsim, f)
    print simpath


#********************************************
# Histogram
#--------------------------------------------
vnameobs  = 'relv'
vnamesim  = 'vortlw'
BBox = [[5,100],[60,180]]
#BBox = [[0,100],[60,180]]
#BBox = [[25,100],[50,180]]
#BBox = [[50,100],[60,180]]
[[lllat,lllon],[urlat,urlon]] = BBox
a1bnd = (np.arange(0, 150+0.01, 1) -0.05) *1e-5
a1cent= 0.5* (a1bnd[:-1] + a1bnd[1:])
#-- Obs ----------
with open(obspath, 'rb') as f:
    dobs = pickle.load(f)
a1varobs = dobs[vnameobs]
a1latobs = dobs['Lat']
a1lonobs = dobs['Lon']

a1mask1 = ma.masked_outside(a1latobs, lllat, urlat).mask
a1mask2 = ma.masked_outside(a1lonobs, lllon, urlon).mask
a1mask  = a1mask1 + a1mask2
a1vartmp = ma.masked_where(a1mask, a1varobs).compressed()

a1histobs, a1bnd = np.histogram(a1vartmp, bins=a1bnd)
a1histobs = a1histobs / float(a1histobs.sum())
#-- Sim ----------
a2histsim = []
for ens in lens:
    simpath = simdir + '/%s.%s.%03d.pickle'%(vnamesim,scen,ens)
    with open(simpath, 'rb') as f:
        dsim = pickle.load(f) 

    a1latsim = dsim['Lat']
    a1lonsim = dsim['Lon']
    a1var    = dsim[vnamesim]
    #-- mask ---
    a1mask1 = ma.masked_outside(a1latsim, lllat, urlat).mask
    a1mask2 = ma.masked_outside(a1lonsim, lllon, urlon).mask
    a1mask  = a1mask1 + a1mask2

    a1vartmp = ma.masked_where(a1mask, a1var).compressed()

    a1histsim, a1bnd = np.histogram(a1vartmp, bins=a1bnd)
    a2histsim.append(a1histsim)

a2histsim = np.array(a2histsim)
a2histsim = a2histsim / a2histsim.sum(axis=1).reshape(-1,1).astype('float32')

#-- Figure ---------------
fig = plt.figure(figsize=(6,5))
ax  = fig.add_axes([0.1,0.1, 0.8, 0.8])

for iens, ens in enumerate(lens):
    print 'plot d4PDF',ens
    a1histsim = a2histsim[iens]
    ax.plot(a1cent, a1histsim, '-', linewidth=1, color='gray')

a1histmedian = np.median(a2histsim, axis=0)
ax.plot(a1cent, a1histmedian, '-', linewidth=1, color='k')


ax.plot(a1cent, a1histobs, '-',linewidth=3, color='k')

stitle = slabel + ' [Lat:%d-%d Lon %d-%d]'%(lllat,urlat,lllon,urlon)
stitle = stitle + '\nUpscale + findmax'
plt.title(stitle)
figpath = '/home/utsumi/temp/bams2020/pdf.%s.lat-%d-%d.lon-%d-%d.png'%(slabel, lllat,urlat, lllon, urlon)
plt.savefig(figpath)
#plt.show()

#--- percentiles (OBS) -------------------------
lpercent = np.arange(0,100+0.1, 1)
aptileobs = np.percentile(a1varobs, lpercent)
lout = []
for i,percent in enumerate(lpercent):
    lout.append([percent, aptileobs[i]])
sout = util.list2csv(lout)
csvpath = '/home/utsumi/temp/bams2020/percentile.%s.tc.obs.lat-%d-%d.lon-%d-%d.csv'%(vnameobs, lllat,urlat, lllon, urlon)
f=open(csvpath, 'w'); f.write(sout); f.close()
print csvpath
# %%
