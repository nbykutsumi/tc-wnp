# %%
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
%matplotlib inline
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import numpy as np
import scipy
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from bisect import bisect_left
from importlib import import_module
from numpy import ma
import APHRODITE

#calcflag = True
calcflag = False
figflag = True
iY = 1951
eY = 2015
lY = range(iY,eY+1)

miss_out = -9999.
hostname = socket.gethostname()
if hostname =='shui':
    dbbaseDir = "/tank/utsumi/data/APHRO/APHRO_V1101"
elif hostname=='well':
    dbbaseDir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101"

region= "MA"
res   = "050"
aph = APHRODITE.v1101(region=region, res=res, dbbaseDir=dbbaseDir)
a1lat = aph.Lat
a1lon = aph.Lon
a1latbnd = aph.LatBnd
a1lonbnd = aph.LonBnd
ny = len(a1lat)  # 140
nx = len(a1lon)  # 180


for Year in lY:
    if calcflag != True: continue

    a3prec, a3site = aph.load_year(Year)
    ndays = a3prec.shape[0]

    a2days_valid = ma.masked_less(a3prec,0).count(axis=0).astype("int16")
    a2days_site  = ma.masked_less_equal(a3site,0).count(axis=0).astype("int16")
    #--- Save -----
    datdir = "/home/utsumi/temp/bams2020/aphro"
    util.mk_dir(datdir)
    validpath = datdir + "/days.valid.%04d.npy"%(Year)
    sitepath  = datdir + "/days.site.%04d.npy"%(Year)
    np.save(validpath, a2days_valid)
    np.save(sitepath,  a2days_site)
    print(a2days_valid.max())
    print(validpath)
 
###***** Figure *********
#lieY = [[1951,1960],[1961,1965],[1966,1970],[1971,1975],[1976,1980],[1981,1985],[1986,1990],[1991,1995],[1996,2000],[2001,2005],[2006,2010],[2011,2015]]


lieY = [[1951,1959],[1960,1969],[1970,1979],[1980,1989],[1990,1999],[2000,2009],[2010,2015]]

a2mask = np.array([np.load(datdir + "/days.valid.%04d.npy"%(Year)).astype("int32") for Year in range(1951,2015+1)]).sum(axis=0)

for (iytmp,eytmp) in lieY:
    lYearTmp = range(iytmp,eytmp+1)
    if figflag != True: continue

    datdir = "/home/utsumi/temp/bams2020/aphro"

    a2valid = np.array([np.load(datdir + "/days.valid.%04d.npy"%(Year)).astype("int32") for Year in lYearTmp]).sum(axis=0)
    a2site  = np.array([np.load(datdir + "/days.site.%04d.npy"%(Year)).astype("int32") for Year in lYearTmp]).sum(axis=0)

    ndays = (datetime(eytmp,12,31) - datetime(iytmp,1,1)).days + 1

    a2valid = a2valid.astype("float32") / ndays
    a2site  = a2site.astype("float32") / ndays

    a2valid = ma.masked_where(a2mask==0, a2valid)
    a2site  = ma.masked_where(a2mask==0, a2site)

    for idat in range(2):
        print(idat)
        if idat==0:
            a2fig = a2valid
        else:
            a2fig = a2site
        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
        axmap.set(facecolor="gray")

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
        cmbnd = list(np.arange(0,1+0.01, 0.1))
        cmlabels = map(float, ["%.1f"%x for x in cmbnd])



        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
        #--extent and coast lines ------
        regionfig = "WNP"
        dbboxfig = {"WNP":[[10,105],[45,150]]}
        [[lllatfig,lllonfig],[urlatfig,urlonfig]] = dbboxfig[regionfig]
        axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
        axmap.coastlines()

        #-- title, figure name -----
        #figdir  = "/home/utsumi/temp/bams2020/fig/map-prec"
        #figpath = figdir + "/map.prec.rp-%03d.APHRO.%s.%04d-%04d.png"%(rp, regionfig, iY, eY)

        stitle = '%s %% %04d-%04d'%({0:"valid",1:"site"}[idat],iytmp,eytmp)

        axmap.set_title(stitle)

        #-- pcolormesh -------------------
        X,Y = np.meshgrid(a1lonbnd, a1latbnd)
        im = axmap.pcolormesh(X,Y, a2fig, cmap=cmap, norm=norm)
        plt.colorbar(im)
        plt.show()

        #plt.savefig(figpath)
        #print(figpath)


# %%
