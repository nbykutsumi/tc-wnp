# %%
import matplotlib
matplotlib.use('Agg')
#----------------------------------
import sys, os, pickle, calendar
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
#--------------------------------------
calcbst= True
#calcbst= False
figbst = True
#figbst = False

calcobj= True
#calcobj= False
#figobj = True
figobj = False

figmean = True
#figmean = False

ctype = 'TC'
#ctype = 'ALL'

lrvort = [[90, 4.0e-4, 7.6e-4],[50, 1.7e-4, 2.0e-4], [0, 4.7e-5, 4.8e-5]] # [ptile(TC), JRA55, d4PDF]
#lrvort = [[90, 4.0e-4, 7.6e-4],[50, 1.7e-4, 2.0e-4]] # [ptile(TC), JRA55, d4PDF]
#lrvort = [[0,4.7e-5, 4.8e-5]] # [ptile(TC), JRA55, d4PDF]
#lrvort = [[50, 1.7e-4, 2.0e-4]] # [ptile(TC), JRA55, d4PDF]

iDTime = datetime(2000,1,1,0)
#eDTime = datetime(2000,1,31,18)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2010,9,1,0)
#eDTime = datetime(2010,9,31,18)

dDTime = timedelta(hours=6)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]

cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = None
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
scen    = 'HPB' # run={expr}-{scen}-{ens}
lens    = range(1,9+1)
#lens    = [1]
res     = "320x640"
noleap  = False


detectName = 'wsd_d4pdf_20200428'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))


wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outDir = '/home/utsumi/mnt/lab_tank/utsumi/d4PDF.EAsia/tune' 
util.mk_dir(outDir)
figdir  = '/home/utsumi/temp/bams2020'
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
[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]
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
#iDTimeBst = datetime(2000,1,1,0)
#eDTimeBst = datetime(2010,12,31)
iYM = [2000,1]
eYM = [2010,12]
#eYM = [2000,12]

lYMbst = util.ret_lYM(iYM, eYM)
for [percent, thrvort, _] in lrvort:
    if calcbst is True:
        dbst   = {}
        bst    = IBTrACS.IBTrACS()

        a2count = np.zeros([ny,nx],int32) 
        for (Year,Mon) in lYMbst:
            print 'Bst',percent,Year,Mon

            eDay   = calendar.monthrange(Year,Mon)[1]
            iDTimetmp = datetime(Year,Mon,1,0)
            eDTimetmp = datetime(Year,Mon,eDay,18)

            #-- Read var ---
            obsvarpath = '/home/utsumi/temp/bams2020/tc-obs-var/relv/%04d.%02d.pickle'%(Year,Mon)
            with open(obsvarpath,'r') as f:
                dvar = pickle.load(f)
            lobsvar = []
            for ldat in dvar.values():
                lobsvar = lobsvar + ldat

            #-- Read besttrack --
            dlonlat = bst.ret_dlonlat(iDTimetmp,eDTimetmp)
            llonlat = []
            for lonlat in dlonlat.values():
                llonlat = llonlat + lonlat


            #-- Map --------
            for i,(lon,lat) in enumerate(llonlat):
                varval = lobsvar[i]
                if varval < thrvort: continue

                ix = int((lon-lonbnd0)/dgrid)
                iy = int((lat-latbnd0)/dgrid)
                if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
                a2count[iy,ix] = a2count[iy,ix] + 1


        iYear,iMon = iYM
        eYear,eMon = eYM
        iDTimeBst = datetime(iYear,iMon,1, 0)
        eDTimeBst = datetime(eYear,eMon,calendar.monthrange(eYear,eMon)[1], 18)
        lDTimeBst = util.ret_lDTime(iDTimeBst, eDTimeBst, timedelta(hours=6))

        a2freq = a2count.astype('float32')/len(lDTimeBst)
        dbst['a2count'] = a2count
        dbst['a2freq' ] = a2freq
        dbst['iDTime']= iDTimeBst
        dbst['eDTime']= eDTimeBst
        dbst['a1latbnd'] = a1latbnd
        dbst['a1lonbnd'] = a1lonbnd
        dbst['dgrid']    = dgrid
        dbst['thrvort']  = thrvort
        dbst['percentile'] = percent

    #-- Save --------
    bstPath= outDir + '/freq.tc.bst.p%02d.pickle'%(percent)
    if calcbst is True:
        with open(bstPath,'wb') as f:
            pickle.dump(dbst, f)
        print bstPath

    #-- Figure (best track)
    if figbst is True:
        with open(bstPath, 'rb') as f:
            dbst = pickle.load(f)

        a2fig = dbst['a2freq']
        a2fig = ma.masked_less(a2fig,0)*4*365  # times/year
        dpara = {}
        dpara['title'] = 'times/year (Best track)' + '\nrvort-p%02d'%(percent)

        figdir  = '/home/utsumi/temp/bams2020'
        dpara['figpath'] = figdir + '/map.freq.tc.bst.p%02d.png'%(percent)
        dpara['cmbnd']   = cmbnd

        draw_map(a2fig, dpara)
#************************************
# d4PDF (Objective detection)
#************************************
#lexrvort= 3.7*1.0e-5 * np.array([1, 1.3])
#lthwcore= [0.2]

#lexrvort= np.array([2.0, 2.5, 3.0])*1.0e-5
#lthwcore= [0.2]

#lthsst  = [25,27,29]
lthsst  = [27]
#lthsst  = [15,25]
#lthsst  = [-273.15]
#lexrvort= np.array([3.7])*1.0e-5
#lexrvort= np.array([2.5,3.0])*1.0e-5
#lexrvort= np.array([2.0, 2.5])*1.0e-5
#lexrvort= np.array([-9999])*1.0e-5
#lexrvort= np.array([6.0])*1.0e-5
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([3,4,5,6])*1.0e-5
lthwcore= [0.3]
#lthwcore= [0.0, 0.2]
#lthwcore= [0, 0.5]
#lthwcore= [-99999]
#lthwcore= [-99]
#lthdura = [48,60,72]
lthdura = [48]
#lthdura = [6,12,24,36]
#lthdura = [-99999,-6,12,24,36]

lkey = [[thsst,exrvort,thwcore,thdura]
        for thsst   in lthsst
        for exrvort in lexrvort
        for thwcore in lthwcore
        for thdura in lthdura
        ]

for (percent,_,tcrvort) in lrvort:
    for (thsst,exrvort,thwcore,thdura) in lkey:
        const  = ConstCyclone.Const(prj=prj, model=model)
        const['Lat'] = d4PDF.Lat()
        const['Lon'] = d4PDF.Lon()
    
        const['thsst']   = thsst + 273.15   # K
        const['exrvort'] = exrvort
        const['tcrvort'] = tcrvort 
        const['thwcore'] = thwcore
        const['thdura']  = thdura
        exrvortout = exrvort*1.0e+5
        tcrvortout = tcrvort*1.0e+5
        slabel = '%s-%s-%s-sst-%d.ex-%.2f.tc-p%02d-%.2f.wc-%.1f-du-%02d'%(ctype, expr, scen, thsst, exrvortout, percent, tcrvortout, thwcore, thdura)
        #slabel = 'wceach-%s-%s-%s-sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(ctype, expr, scen, thsst, exrvortout, tcrvortout, thwcore, thdura)
    
        for ens in lens:
            if calcobj is True:
                ltclonlat = []
                print 'ens=',ens
                #-------------------
                wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    
                cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)
                dexcxy, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
    
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
                objPath= outDir + '/freq.tc.obj.%s.en-%03d.pickle'%(slabel, ens)
                if calcobj is True:
                    with open(objPath,'wb') as f:
                        pickle.dump(dobj, f)
                    print objPath
    
            #------------------------
            # Figure for each ensemble
            objPath= outDir + '/freq.tc.obj.%s.en-%03d.pickle'%(slabel, ens)
            if figobj is True:
                with open(objPath, 'rb') as f:
                    dobj = pickle.load(f)
    
                a2fig = dobj['a2freq']
                a2fig = ma.masked_less(a2fig,0)*4*365  # times/10-year
    
                dpara = {}
                dpara['title'] = 'time/year (d4PDF)' + '\n' + '%s-%s sst:%d ex:%.2f tc:p%02d-%.2f wc:%.1f ens:%03d'%(expr, scen, thsst, exrvortout, percent, tcrvortout, thwcore, ens)
                #dpara['figpath'] = figdir + '/map.freq.tc.obj.png'
                dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.en-%03d.png'%(slabel, ens)
                dpara['cmbnd'] = cmbnd
                draw_map(a2fig, dpara)
    
        #***********************
        if figmean is True:
            a2count = np.zeros([ny,nx],int32) 
            nstep   = 0
            for ens in lens:
                objPath= outDir + '/freq.tc.obj.%s.en-%03d.pickle'%(slabel, ens)
                with open(objPath, 'rb') as f:
                    dobj = pickle.load(f)
    
                a2countTmp = dobj['a2count']
                nstepTmp   = dobj['nstep']
                a2count    = a2count + a2countTmp
                nstep      = nstep + nstepTmp
    
            a2fig = a2count.astype('float32') / nstep
            a2fig = ma.masked_less(a2fig,0)*4*365  # times/10-year
    
            dpara = {}
            dpara['title'] = 'times/year (d4PDF) ' + '\n' + '%s-%s sst:%d ex:%.2f tc:p%02d-%.2f wc:%.1f ens-mean'%(expr, scen, thsst, exrvortout, percent, tcrvortout, thwcore)
            dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.ave.png'%(slabel)
            dpara['cmbnd'] = cmbnd
            draw_map(a2fig, dpara)

# %%
