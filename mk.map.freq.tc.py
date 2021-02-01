# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
import socket
#import Cyclone
#--------------------------------------
#calcbst= True
calcbst= False
figbst = True
#figbst = False

#calcobj= True
calcobj= False
#figobj = True
figobj = False

figmean = True
#figmean = False

ctype = 'TC'
#ctype = 'ALL'

iY, eY = 1980,2010
#iY, eY = 1990,2010
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])

#cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = np.arange(0,50,4)
#cmbnd = [0,1,5] + range(10,50,5)
#cmbnd = np.array(cmbnd)
#cmbnd = np.arange(0,50,4)
cmbnd = np.arange(0,50,4)
print(cmbnd)
#cmbnd = None
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB_NAT' # run={expr}-{scen}-{ens}
lens    = range(1,20+1)
#lens    = range(21,50+1)
#lens    = range(1,50+1)
#lens    = [20]
#lens    = range(3,9+1)
res     = "320x640"
noleap  = False


detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))  # test
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))


wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outbaseDir = '/home/utsumi/temp/bams2020/map-freq'
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-freq'
util.mk_dir(figdir)
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999
dgrid = 5
a1latbnd  = np.arange(-90,90+0.01, dgrid)
a1lonbnd  = np.arange(0,360+0.01, dgrid)
ny = len(a1latbnd) - 1
nx = len(a1lonbnd) - 1
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]
[[lllat,lllon],[urlat,urlon]] = [[10,105],[45,150]]
a1latbndfig = np.arange(lllat, urlat+0.01, dgrid)
a1lonbndfig = np.arange(lllon, urlon+0.01, dgrid)
#nyfig = a1latbnd.shape[0] -1
#nxfig = a1lonbnd.shape[0] -1

#**************************
def draw_map(a2dat, dpara):
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())

    gl        = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,901, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)

    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines(color="k")

    ##-- Make new colormap (white at the lower end) --
    #upper = matplotlib.cm.jet(np.arange(256))
    #lower = np.ones((int(256/4),4))
    #for i in range(3):
    #  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
    #mycm = np.vstack(( lower, upper ))
    #mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

    mycm = "gist_stern_r"

    #-- color boundaries norm --------
    #mycm = 'gist_stern_r'
    cmbnd = dpara['cmbnd']
    if cmbnd is not None:
        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]

        cmap = matplotlib.colors.ListedColormap(cmaplist)

        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

    else:
        cmap = mycm
        norm = None


    #-- pcolormesh --------------
    vmin, vmax = 0, None
    X,Y = np.meshgrid(a1lonbnd,a1latbnd)
    im  = plt.pcolormesh(X,Y,a2dat, cmap=cmap, vmin=vmin,vmax=vmax, norm=norm)
    #-- Colorbar ------
    cax = figmap.add_axes([0.82,0.2,0.05,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)

    if cmbnd is not None:
        extend = dpara['extend']
        cbar.set_ticks(cbar.ax.get_yticks())
        if extend=='both':
            cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:-1] + [""])

        if extend=='min':
            cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:])

        if extend=='max':
            cbar.set_ticklabels(list(cbar.ax.get_yticklabels())[:-1]+ [""])

        else:
            pass
    #-- coastline ---------------

    #-- Tiele -----
    stitle = dpara['title']
    axmap.set_title(stitle)
    #-- Save ------
    figpath = dpara['figpath']
    plt.savefig(figpath)
    plt.show()
    print(figpath)



#**************************
#Read best track
#---------------------------
#iYbst,eYbst = 2000, 2010
#iYbst,eYbst = 1980, 2018
iYbst,eYbst = 1980, 2010
#iYbst,eYbst = 1990, 2010
#iYbst,eYbst = 2018,2018


lYMbst = util.ret_lYM([iYbst,1],[eYbst,12])

for (Year,Mon) in lYMbst:
    YM = (Year,Mon)
    if calcbst is not True: continue
    print((Year,Mon))

    lDTimeBst = util.ret_lDTime_fromYM(YM, YM, timedelta(hours=6), 0)
    iDTimeBst = lDTimeBst[0]
    eDTimeBst = lDTimeBst[-1]
 
    dbst   = {}
    bst    = IBTrACS.IBTrACS()
    dlonlat = bst.ret_dlonlat(iDTimeBst,eDTimeBst)
    llonlat = []
    for lonlat in list(dlonlat.values()):
        llonlat = llonlat + lonlat
    #-- Map --------
    a2count = np.zeros([ny,nx],int32) 
    for (lon,lat) in llonlat:
        ix = int((lon-lonbnd0)/dgrid)
        iy = int((lat-latbnd0)/dgrid)
        if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
        a2count[iy,ix] = a2count[iy,ix] + 1

    a2freq = a2count.astype('float32')/len(lDTimeBst)

    dbst['a2count'] = a2count
    dbst['a2freq' ] = a2freq
    dbst['iDTime']= iDTimeBst
    dbst['eDTime']= eDTimeBst
    dbst['nstep' ]= len(lDTimeBst)
    dbst['a1latbnd'] = a1latbnd
    dbst['a1lonbnd'] = a1lonbnd
    dbst['dgrid']    = dgrid

    #-- Save --------
    bstdir = outbaseDir + '/bst'
    util.mk_dir(bstdir)

    for vname in list(dbst.keys()):
        bstPath= bstdir + '/%s.tc.bst.%04d.%02d.npy'%(vname,Year,Mon)
        np.save(bstPath, dbst[vname] ) 
        print(bstPath)

#-- Figure (best track)
if figbst is True:
    a2count = np.zeros([ny,nx], 'int32')
    nstep   = 0
    for (Year,Mon) in lYMbst:
        bstdir = outbaseDir + '/bst'
    
        a2countTmp = np.load(bstdir + '/a2count.tc.bst.%04d.%02d.npy'%(Year,Mon))
        nstepTmp   = np.load(bstdir + '/nstep.tc.bst.%04d.%02d.npy'%(Year,Mon))

        a2count = a2count + a2countTmp
        nstep   = nstep   + nstepTmp
    
    a2fig = a2count.astype('float32') / nstep
    a2fig = ma.masked_less(a2fig,0)*4*365  # times/year
    dpara = {}
    dpara['title'] = 'count/10-year (Best track) %04d-%04d'%(iYbst,eYbst)
    dpara['figpath'] = figdir + '/map.freq.tc.bst.%04d-%04d.png'%(iYbst,eYbst)
    dpara['cmbnd']   = cmbnd
    dpara['extend']  = "max"
    draw_map(a2fig, dpara)

#************************************
# d4PDF (Objective detection)
#************************************
#lthsst  = [27,28]
#lthsst  = [27,27.5,28]
lthsst  = [27]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
#lthwcore= [-1,0,1]
lthdura = [36]
#lthwind = [10,13,15]
lthwind = [12]
#lthwdif = [-9999]
lthwdif = [-9999]

lkey = [[thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif]
        for thsst   in lthsst
        for exrvort in lexrvort
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thdura in lthdura
        for thwind in lthwind
        for thwdif in lthwdif
        ]


for (thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif) in lkey:
    const  = ConstCyclone.Const(prj=prj, model=model)
    const['Lat'] = d4PDF.Lat()
    const['Lon'] = d4PDF.Lon()

    const['thsst']   = thsst + 273.15   # K
    const['exrvort'] = exrvort
    const['tcrvort'] = tcrvort 
    const['thwcore'] = thwcore
    const['thdura']  = thdura
    const['thwind']  = thwind
    const['thwdif']  = thwdif

    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    for ens in lens:
        for Year in lYear:
            if calcobj is not True: continue
            for Mon in lMon:

                eday   = calendar.monthrange(Year,Mon)[1]
                iDTime = datetime(Year,Mon,1,0)
                eDTime = datetime(Year,Mon,eday,18)
                lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=6))

                ltclonlat = []
                print(('ens=',ens))
                #-------------------
                wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    
                cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)
                dexcxy, dtcxy  = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname='vortlw')

                ltcxy = []
                for ltmp in list(dtcxy.values()):
                    ltcxy = ltcxy + ltmp

                #-------------
                if ctype == 'ALL': 
                    for ltmp in list(dexcxy.values()):
                        ltcxy = ltcxy + ltmp
                #-------------
                for tcxy in ltcxy:
                    x,y,var = tcxy
                    lon,lat = a1lonin[x],a1latin[y]
                    ltclonlat.append([lon,lat,var])


                #-- Map --------
                a2count = np.zeros([ny,nx],int32) 
                for (lon,lat,_) in ltclonlat:
                    ix = int((lon-lonbnd0)/dgrid)
                    iy = int((lat-latbnd0)/dgrid)
                    if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
                    a2count[iy,ix] = a2count[iy,ix] + 1

                a2freq = a2count.astype('float32')/len(lDTime)

                dobj   = {}
                dobj['a2count'] = a2count
                dobj['a2freq' ] = a2freq
                dobj['iDTime']= iDTime
                dobj['eDTime']= eDTime
                dobj['nstep'] = len(lDTime)
                dobj['a1latbnd'] = a1latbnd
                dobj['a1lonbnd'] = a1lonbnd
                dobj['dgrid']    = dgrid
                dobj['const'] = const
                dobj['tcrvort'] = tcrvort

                #-- Save --------
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                util.mk_dir(outDir)

                for vname in dobj.keys():
                    objPath= outDir + '/%s.tc.obj.%04d.%02d.npy'%(vname,Year,Mon)
                    np.save(objPath, dobj[vname])
                    print(objPath)

            #-- Make annual data -----
            a2count = np.zeros([ny,nx],int32) 
            nstep   = 0
            outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
            for Mon in lMon:
                a2countTmp = np.load(outDir + '/a2count.tc.obj.%04d.%02d.npy'%(Year,Mon))
                nstepTmp   = np.load(outDir + '/nstep.tc.obj.%04d.%02d.npy'%(Year,Mon))

                a2count = a2count + a2countTmp
                nstep   = nstep + nstepTmp

            a2freq = a2count.astype("float32") / nstep

            dobj   = {}
            dobj['a2count'] = a2count
            dobj['a2freq' ] = a2freq
            dobj['iDTime']= datetime(Year,1,1,0)
            dobj['eDTime']= datetime(Year,12,31,18)
            dobj['nstep'] = nstep
            dobj['a1latbnd'] = a1latbnd
            dobj['a1lonbnd'] = a1lonbnd
            dobj['dgrid']    = dgrid
            dobj['const'] = const
            dobj['tcrvort'] = tcrvort

            for vname in dobj.keys():
                objPath= outDir + '/%s.tc.obj.%04d.npy'%(vname,Year)
                np.save(objPath, dobj[vname])
                print(objPath)


        #------------------------
        # Figure for each ensemble
        a2count = np.zeros([ny,nx],'int32')
        nstep   = 0

        if figobj is True:
            for Year in lYear:
    
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                a2countTmp = np.load(outDir + '/a2count.tc.obj.%04d.npy'%(Year,Mon))
                nstepTmp   = np.load(outDir + '/nstep.tc.obj.%04d.npy'%(Year,Mon))
 
                a2count = a2count + a2countTmp
                nstep   = nstep   + nstepTmp
    
    
            a2fig = a2count.astype('float32')/nstep
            a2fig = ma.masked_less(a2fig,0)*4*365  # times/10-year
    
            dpara = {}
            dpara['title'] = 'count/10-year (d4PDF)' + '%04d-%04d'%(iY,eY) + '\n' + '%s-%s sst:%d ex:%.2f tc:%.2f \n wc:%.1f wind:%d wdif:%d ens:%03d'%(expr, scen, thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, ens)
            #dpara['figpath'] = figdir + '/map.freq.tc.obj.png'
            dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.%04d-%04d.en-%03d.png'%(slabel, iY, eY, ens)
            dpara['cmbnd'] = cmbnd
            #draw_map(a2fig, dpara)

    #***********************
    # Figure: ensemble mean
    #***********************
    a2count = np.zeros([ny,nx],'int32')
    nstep   = 0

    if figmean is True:
        for ens in lens:
            for Year in lYear:
                print(ens,Year) 
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                print(outDir + '/a2count.tc.obj.%04d.npy'%(Year))
                a2countTmp = np.load(outDir + '/a2count.tc.obj.%04d.npy'%(Year))
                nstepTmp   = np.load(outDir + '/nstep.tc.obj.%04d.npy'%(Year))
    
                a2count = a2count + a2countTmp
                nstep   = nstep   + nstepTmp

        a2fig = a2count.astype('float32') / nstep
        a2fig = ma.masked_less(a2fig,0)*4*365  # times/10-year
    
        dpara = {}
        dpara['title'] = 'Prob. of existence (d4PDF) %04d-%04d'%(iY, eY) + '\n' + '%s-%s sst:%d ex:%.2f tc:%.2f \n wc:%.1f wind:%d wdif:%d ens-mean'%(expr, scen, thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif)
        dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.%04d-%04d.ave.png'%(slabel, iY, eY)
        dpara['cmbnd'] = cmbnd
        dpara['extend']= "max"
        draw_map(a2fig, dpara)

# %%
