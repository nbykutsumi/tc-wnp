# %%
import matplotlib
matplotlib.use('Agg')
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
import calendar
#--------------------------------------
#calcbst= True
calcbst= False
figbst = True
#figbst = False

#calcobj= True
calcobj= False
figobj = True
#figobj = False

figmean = True
#figmean = False

ctype = 'TC'
#ctype = 'ALL'

iY, eY = 2000,2010
#iY, eY = 2001,2010
lYM = util.ret_lYM([iY,1],[eY,12])

#cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = np.arange(0,50,4)
#cmbnd = [0,1,5] + range(10,50,5)
#cmbnd = np.array(cmbnd)
cmbnd = np.arange(0,50,4)
print cmbnd
#cmbnd = None
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB_NAT' # run={expr}-{scen}-{ens}
#lens    = range(1,9+1)
lens    = range(1,2+1)
#lens    = range(3,9+1)
res     = "320x640"
noleap  = False


detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
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
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]
a1latbndfig = np.arange(lllat, urlat+0.01, dgrid)
a1lonbndfig = np.arange(lllon, urlon+0.01, dgrid)
#nyfig = a1latbnd.shape[0] -1
#nxfig = a1lonbnd.shape[0] -1

#**************************
def draw_map(a2dat, dpara):
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8])
    M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

    #-- color boundaries norm --------
    mycm = 'gist_stern_r'
    bounds = dpara['cmbnd']
    if bounds is not None:
        cmap   = plt.cm.get_cmap(mycm, len(bounds)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]

        cmap = matplotlib.colors.ListedColormap(cmaplist)

        norm = matplotlib.colors.BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

    else:
        cmap = mycm
        norm = None

    #-- pcolormesh --------------
    vmin, vmax = 0, None
    X,Y = np.meshgrid(a1lonbnd,a1latbnd)
    im  = M.pcolormesh(X,Y,a2dat, cmap=cmap, vmin=vmin,vmax=vmax, norm=norm)
    #-- Colorbar ------
    cax = figmap.add_axes([0.82,0.2,0.05,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)

    #if bounds is not None:
    #    extend = dpara['extend']
    #    if extend=='both':
    #        cbar.ax.set_yticklabels([''] + list(bounds[1:-1]) + [''])
    #    if extend=='min':
    #        cbar.ax.set_yticklabels([''] + list(bounds[1:]))
    #    if extend=='max':
    #        cbar.ax.set_yticklabels(list(bounds[:-1]) + [''])

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
#iYbst,eYbst = 2000, 2010
iYbst,eYbst = 1980, 2018
#iYbst,eYbst = 2018,2018


lYMbst = util.ret_lYM([iYbst,1],[eYbst,12])

for (Year,Mon) in lYMbst:
    YM = (Year,Mon)
    if calcbst is not True: continue
    print Year,Mon

    lDTimeBst = util.ret_lDTime_fromYM(YM, YM, timedelta(hours=6), 0)
    iDTimeBst = lDTimeBst[0]
    eDTimeBst = lDTimeBst[-1]
 
    dbst   = {}
    bst    = IBTrACS.IBTrACS()
    dlonlat = bst.ret_dlonlat(iDTimeBst,eDTimeBst)
    llonlat = []
    for lonlat in dlonlat.values():
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

    bstPath= bstdir + '/freq.tc.bst.%04d.%02d.pickle'%(Year,Mon)
    with open(bstPath,'wb') as f:
        pickle.dump(dbst, f)
    print bstPath


#-- Figure (best track)
if figbst is True:
    a2count = np.zeros([ny,nx], 'int32')
    nstep   = 0
    for (Year,Mon) in lYMbst:
        bstdir = outbaseDir + '/bst'
        bstPath= bstdir + '/freq.tc.bst.%04d.%02d.pickle'%(Year,Mon)
    
        with open(bstPath, 'rb') as f:
            dbst = pickle.load(f)
    
        a2count = a2count + dbst['a2count']
        nstep   = nstep   + dbst['nstep']
    
    a2fig = a2count.astype('float32') / nstep
    a2fig = ma.masked_less(a2fig,0)*4*365  # times/year
    dpara = {}
    dpara['title'] = 'count/10-year (Best track)'
    dpara['figpath'] = figdir + '/map.freq.tc.bst.%04d-%04d.png'%(iYbst,eYbst)
    dpara['cmbnd']   = cmbnd
    
    draw_map(a2fig, dpara)

#************************************
# d4PDF (Objective detection)
#************************************
#lthsst  = [27,28]
#lthsst  = [27,27.5,28]
lthsst  = [27,27.5, 28]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
#lthwcore= [-1,0,1]
lthdura = [36]
#lthwind = [10,13,15]
lthwind = [10,11,12,13]
#lthwdif = [-9999]
lthwdif = [-9999,-30,-20,-10]

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
        for (Year,Mon) in lYM:
            if calcobj is not True: continue

            eday   = calendar.monthrange(Year,Mon)[1]
            iDTime = datetime(Year,Mon,1,0)
            eDTime = datetime(Year,Mon,eday,18)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=6))

            ltclonlat = []
            print 'ens=',ens
            #-------------------
            wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
 
            cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)
            dexcxy, dtcxy  = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname='vortlw')

            ltcxy = []
            for ltmp in dtcxy.values():
                ltcxy = ltcxy + ltmp

            #-------------
            if ctype is 'ALL': 
                for ltmp in dexcxy.values():
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
            objPath= outDir + '/freq.tc.obj.en-%03d.%04d.%02d.pickle'%(ens,Year,Mon)
            if calcobj is True:
                with open(objPath,'wb') as f:
                    pickle.dump(dobj, f)
                print objPath

        #------------------------
        # Figure for each ensemble
        a2count = np.zeros([ny,nx],'int32')
        nstep   = 0

        if figobj is True:
            for (Year,Mon) in lYM:
    
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                objPath= outDir + '/freq.tc.obj.en-%03d.%04d.%02d.pickle'%(ens,Year,Mon)
    
    
                with open(objPath, 'rb') as f:
                    dobj = pickle.load(f)
    
                a2count = a2count + dobj['a2count']
                nstep   = nstep   + dobj['nstep']
    
    
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
            for (Year,Mon) in lYM:
    
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                objPath= outDir + '/freq.tc.obj.en-%03d.%04d.%02d.pickle'%(ens,Year,Mon)
    
                with open(objPath, 'rb') as f:
                    dobj = pickle.load(f)
    
                a2count = a2count + dobj['a2count']
                nstep   = nstep   + dobj['nstep']
    
        a2fig = a2count.astype('float32') / nstep
        a2fig = ma.masked_less(a2fig,0)*4*365  # times/10-year
    
        dpara = {}
        dpara['title'] = 'Prob. of existence (d4PDF) %04d-%04d'%(iY, eY) + '\n' + '%s-%s sst:%d ex:%.2f tc:%.2f \n wc:%.1f wind:%d wdif:%d ens-mean'%(expr, scen, thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif)
        dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.%04d-%04d.ave.png'%(slabel, iY, eY)
        dpara['cmbnd'] = cmbnd
        draw_map(a2fig, dpara)

# %%
