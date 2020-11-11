# %%
import matplotlib
matplotlib.use('Agg')
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   detect_fsub import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
#--------------------------------------
#calcbst= True
calcbst= False
#figbst = True
figbst = False

calcobj= True
#calcobj= False
figobj = True
#figobj = False

iDTime = datetime(2006,6,1,0)
eDTime = datetime(2015,12,31,0)
#eDTime = datetime(2008,12,31,0)
dDTime = timedelta(hours=6)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
#-----------------
prj     = "HAPPI"
model   = "MIROC5"
runpart = "C20-ALL"
res     = "128x256"
lens    = range(1,50+1)
#lens    = range(1,1+1)

detectName = 'detect_20200401'
config_func = import_module("%s.config_func"%(detectName))
ConstC      = import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
IO_Master   = import_module("%s.IO_Master"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))

outDir = '/home/utsumi/mnt/lab_tank/utsumi/HAPPI.EAsia/tune' 
util.mk_dir(outDir)
figdir  = '/home/utsumi/temp/happi'
#----------------------------------
iom = IO_Master.IO_Master(prj, model, 'C20-ALL-001', res)
a1latin = iom.Lat
a1lonin = iom.Lon
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999

[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
dgrid = 5
a1latbnd = np.arange(lllat, urlat+0.01, dgrid)
a1lonbnd = np.arange(lllon, urlon+0.01, dgrid)
ny = a1latbnd.shape[0] -1
nx = a1lonbnd.shape[0] -1

#**************************
def draw_map(a2dat, dpara):
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8])
    M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

    #-- pcolormesh --------------
    mycm = 'gist_stern_r'
    vmin, vmax = 0, 0.02
    X,Y = np.meshgrid(a1lonbnd,a1latbnd)
    im  = M.pcolormesh(X,Y,a2dat, cmap=mycm, vmin=vmin,vmax=vmax)
    #-- Colorbar ------
    cax = figmap.add_axes([0.82,0.2,0.05,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)

    #-- coastline ---------------
    print "coastlines"
    M.drawcoastlines(color="k")

    #-- Merdians and Parallels --
    M.drawmeridians(np.arange(-180,180+1,15), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50, yoffset=2)
    M.drawparallels(np.arange(-60,60+1,15), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    #-- Tiele -----
    stitle = dpara['title']
    axmap.set_title(stitle)
    #-- Save ------
    figpath = dpara['figpath']
    plt.savefig(figpath)
    print figpath



#**************************
#Read best track
#---------------------------
iDTimeBst = datetime(1980,1,1,0)
eDTimeBst = datetime(2019,12,31)
lDTimeBst = util.ret_lDTime(iDTimeBst,eDTimeBst,timedelta(hours=6))
if calcbst is True:
    dbst   = {}
    bst    = IBTrACS.IBTrACS()
    dlonlat = bst.ret_dlonlat(iDTimeBst,eDTimeBst)
    llonlat = []
    for lonlat in dlonlat.values():
        llonlat = llonlat + lonlat
    #-- Map --------
    a2count = np.zeros([ny,nx],int32) 
    for (lon,lat) in llonlat:
        ix = int((lon-lllon)/dgrid)
        iy = int((lat-lllat)/dgrid)
        if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
        a2count[iy,ix] = a2count[iy,ix] + 1

    a2freq = a2count.astype('float32')/len(lDTimeBst)

    dbst['a2count'] = a2count
    dbst['a2freq' ] = a2freq
    dbst['iDTime']= iDTimeBst
    dbst['eDTime']= eDTimeBst
    dbst['a1latbnd'] = np.arange(lllat, urlat+0.01, dgrid)
    dbst['a1lonbnd'] = np.arange(lllon, urlon+0.01, dgrid)
    dbst['BBox']     = [[lllat,lllon],[urlat,urlon]]
    dbst['dgrid']    = dgrid

#-- Save --------
bstPath= outDir + '/freq.tc.bst.pickle'
if calcbst is True:
    with open(bstPath,'wb') as f:
        pickle.dump(dbst, f)
    print bstPath

#-- Figure (best track)
if figbst is True:
    with open(bstPath, 'rb') as f:
        dbst = pickle.load(f)

    a2fig = dbst['a2freq']
    dpara = {}
    dpara['title'] = 'Prob. of existence (Best track)'
    figdir  = '/home/utsumi/temp/happi'
    dpara['figpath'] = figdir + '/map.freq.tc.bst.png'

    draw_map(a2fig, dpara)


#************************************
# HAPPI (Objective detection)
#************************************
#lexrvort= 3.7*1.0e-5 * np.array([1, 1.3])
#ltcrvortfactor = [1, 1.3]
#lthwcore= [0.2]

#lexrvort= np.array([2.0, 2.5, 3.0])*1.0e-5
#ltcrvortfactor = [1]
#lthwcore= [0.2]

lexrvort= np.array([4.8])*1.0e-5
ltcrvortfactor = [1]
lthwcore= [0, 0.2]
lthdura = [36]

lkey = [[exrvort,tcrvortfactor,thwcore,thdura]
        for exrvort in lexrvort
        for tcrvortfactor in ltcrvortfactor
        for thwcore in lthwcore
        for thdura in lthdura
        ]

for (exrvort,tcrvortfactor,thwcore,thdura) in lkey:

    if calcobj is True:
        ltclonlat = []
        for ens in lens:
            print 'ens=',ens
            run = runpart + '-%03d'%(ens)
            cfg = config_func.config_func(prj=prj, model=model, run=run)
            const=ConstC.Const(cfg)
            const['exrvort'] = exrvort
            const['tcrvort'] = exrvort * tcrvortfactor
            const['thwcore'] = thwcore

            cy  = Cyclone.Cyclone(cfg=cfg, const=const)

            _, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')

            ltcxy = []
            for ltmp in dtcxy.values():
                ltcxy = ltcxy + ltmp

            for tcxy in ltcxy:
                x,y,var = tcxy
                lon,lat = a1lonin[x],a1latin[y]
                ltclonlat.append([lon,lat,var])

        #-- Map --------
        a2count = np.zeros([ny,nx],int32) 
        for (lon,lat,_) in ltclonlat:
            ix = int((lon-lllon)/dgrid)
            iy = int((lat-lllat)/dgrid)
            if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
            a2count[iy,ix] = a2count[iy,ix] + 1

        a2freq = a2count.astype('float32')/(len(lDTime)*len(lens))

        dobj   = {}
        dobj['a2count'] = a2count
        dobj['a2freq' ] = a2freq
        dobj['iDTime']= iDTime
        dobj['eDTime']= eDTime
        dobj['a1latbnd'] = np.arange(lllat, urlat+0.01, dgrid)
        dobj['a1lonbnd'] = np.arange(lllon, urlon+0.01, dgrid)
        dobj['BBox']     = [[lllat,lllon],[urlat,urlon]]
        dobj['dgrid']    = dgrid
        dobj['cfg'] = cfg
        dobj['const'] = const
        dobj['tcrvortfactor'] = tcrvortfactor
    #-- Save --------
    exrvortout = exrvort*1.0e+5
    slabel = 'ex-%.2f.tc-%.1f.wc-%.1f-du-%02d'%(exrvortout, tcrvortfactor, thwcore, thdura)
    objPath= outDir + '/freq.tc.obj.%s.pickle'%(slabel)
    if calcobj is True:
        with open(objPath,'wb') as f:
            pickle.dump(dobj, f)
        print objPath

    #------------------------
    if figobj is True:
        with open(objPath, 'rb') as f:
            dobj = pickle.load(f)

        a2fig = dobj['a2freq']
        print ''
        print objPath
        print slabel
        print a2fig.mean()
        print dobj['const']
        print ''

        dpara = {}
        dpara['title'] = 'Prob. of existence (HAPPI)' + '\n' + 'ex:%.2f tc:%.1f wc:%.1f'%(exrvortout, tcrvortfactor, thwcore)
        figdir  = '/home/utsumi/temp/happi'
        #dpara['figpath'] = figdir + '/map.freq.tc.obj.png'
        dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.png'%(slabel)
        draw_map(a2fig, dpara)


# %%
