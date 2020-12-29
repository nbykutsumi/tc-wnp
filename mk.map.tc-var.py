# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

#----------------------------------
import sys, os, pickle
#from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
#import Cyclone
#--------------------------------------
#calcobj= True
calcobj= False
figflag = True
#figflag = False

#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,20+1)
#lens    = range(21,50+1)
lens    = range(1,50+1)
#lens    = [20]
#lens    = range(3,9+1)
res     = "320x640"
noleap  = False
#tcvar   = "dslp"
tcvar   = "wmaxlw"

detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))  # test
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))


wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outbaseDir = '/home/utsumi/temp/bams2020/map-tc-%s'%(tcvar)
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-tc-var'
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

    for scen in lscen:
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
                    if tcvar=="dslp":
                        dexcxy, dtcxy_slp = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname="slp")
                        dexcxy, dtcxy_ave = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname="slp_mean_adj")

                        ltcxy = []
                        for ltmp in list(dtcxy_slp.values()):
                            ltcxy = ltcxy + ltmp

                        lvar2= []
                        for ltmp in list(dtcxy_ave.values()):
                            lvar2 = lvar2 + ltmp

                        atcxy = np.array(ltcxy)
                        avar2 = np.array(lvar2)
                        atcxy[:,2] = avar2[:,2] - atcxy[:,2]

                    else:
                        dexcxy, dtcxy = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname=tcvar)

                        ltcxy = []
                        for ltmp in list(dtcxy.values()):
                            ltcxy = ltcxy + ltmp
                        atcxy = np.array(ltcxy)

                    a1x = atcxy[:,0].astype("int32")
                    a1y = atcxy[:,1].astype("int32")              
                    a1v = atcxy[:,2]
                    ltcxy = zip(a1x.tolist(), a1y.tolist(), a1v.tolist())
                    #-------------
                    for tcxy in ltcxy:
                        x,y,var = tcxy
                        lon,lat = a1lonin[x],a1latin[y]
                        ltclonlat.append([lon,lat,var])


                    #-- Map --------
                    a2sum = np.zeros([ny,nx],float64) 
                    a2num = np.zeros([ny,nx],int32) 
                    for (lon,lat,var) in ltclonlat:
                        ix = int((lon-lonbnd0)/dgrid)
                        iy = int((lat-latbnd0)/dgrid)
                        if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
                        a2sum[iy,ix] = a2sum[iy,ix] + var
                        a2num[iy,ix] = a2num[iy,ix] + 1

                    dobj   = {}
                    dobj['a2sum'] = a2sum
                    dobj['a2num' ]= a2num
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
                a2sum = np.zeros([ny,nx],float64) 
                a2num = np.zeros([ny,nx],int32) 
                nstep   = 0
                outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                for Mon in lMon:
                    a2sumTmp = np.load(outDir + '/a2sum.tc.obj.%04d.%02d.npy'%(Year,Mon))
                    a2numTmp = np.load(outDir + '/a2num.tc.obj.%04d.%02d.npy'%(Year,Mon))
                    nstepTmp   = np.load(outDir + '/nstep.tc.obj.%04d.%02d.npy'%(Year,Mon))

                    a2sum = a2sum + a2sumTmp
                    a2num = a2num + a2numTmp
                    nstep   = nstep + nstepTmp


                dobj   = {}
                dobj['a2sum'] = a2sum
                dobj['a2num' ] = a2num
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



###*** Draw map *************
for scen in lscen:
    if figflag !=True: continue

    a3sum = np.empty([len(lens), ny,nx], "float64")
    a3num = np.empty([len(lens), ny,nx], "int32")
    for iens, ens in enumerate(lens):
        print(scen,ens)
        outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)

        a3sum[iens] = np.array([np.load(outDir + "/a2sum.tc.obj.%04d.npy"%(Year)) for Year in lYear]).sum(axis=0)
        a3num[iens] = np.array([np.load(outDir + "/a2num.tc.obj.%04d.npy"%(Year)) for Year in lYear]).sum(axis=0)


    a2sum = a3sum.sum(axis=0)
    a2num = a3num.sum(axis=0)
    a2ave = (ma.masked_where(a2num==0, a2sum) / a2num).filled(-9999.)

    a2fig = a2ave
    if tcvar=="dslp":
        a2fig = a2fig*0.01 # hPa

    a1lat = d4PDF.Lat()
    a1lon = d4PDF.Lon()
    BBox = [[0,100],[50,150]]

    [[lllat,lllon],[urlat,urlon]] = BBox

    fig = plt.figure(figsize=(6,4))
    axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())

    gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks = np.arange(-180,180+1,15)
    yticks = np.arange(-90,90+1,15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
    axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

    #-- Make new colormap (white at the lower end) --
    upper = matplotlib.cm.jet(np.arange(256))
    lower = np.ones((int(256/4),4))
    for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
    mycm = np.vstack(( lower, upper ))
    mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])


    #-- color boundaries norm ------
    if tcvar=="dslp":
        cmbnd = list(np.arange(0,5+0.1,0.25))
        cmlabels = list(np.arange(0,5+0.1,1))
    elif tcvar=="wmaxlw":
        cmbnd = list(np.arange(0,40+0.1,2))
        cmlabels = list(np.arange(0,40+0.1,4))


    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #-- title, figure name -----
    if tcvar=="dslp":
        sunit = "hPa"
    elif tcvar=="wmaxlw":
        sunit = "m/s"

    figpath = figdir + '/map.tc-%s.%s.%04d-%04d.png'%(tcvar,scen,iY,eY)
    stitle = 'TC %s (%s) %s\n'%(tcvar, sunit, scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])

    axmap.set_title(stitle)


    #-- pcolormesh -------------------
    X,Y = np.meshgrid(a1lonbnd, a1latbnd)
    im = axmap.pcolormesh(X,Y, a2fig, cmap=cmap, norm=norm)
    cbar=plt.colorbar(im)
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)

    plt.show()

    plt.savefig(figpath)
    print(figpath)
                               
# %%
