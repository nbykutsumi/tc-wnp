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

figflag = True
#iY,eY = 1960,2015
iY,eY = 1980,2015
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


datdir = "/home/utsumi/temp/bams2020/aphro"
a2site  = np.array([np.load(datdir + "/days.site.%04d.npy"%(Year)).astype("int32") for Year in lY]).sum(axis=0)
ndays = (datetime(eY,12,31) - datetime(iY,1,1)).days + 1

a2site = a2site.astype("float32")/ndays

#-- expand grids with >90% ---
thratio = 0.9
a2wide = np.zeros([ny+2, nx+2], "int16")
a2org  = np.ones([ny,nx],"int16")
a2org[a2site<thratio]=0 
a2wide[0:ny,0:nx]     = a2wide[0:ny,0:nx] + a2org
a2wide[0:ny,1:1+nx]   = a2wide[0:ny,1:1+nx] + a2org
a2wide[0:ny,2:2+nx]   = a2wide[0:ny,2:2+nx] + a2org
a2wide[1:1+ny,0:nx]   = a2wide[1:1+ny,0:nx] + a2org
a2wide[1:1+ny,1:1+nx] = a2wide[1:1+ny,1:1+nx] + a2org
a2wide[1:1+ny,2:2+nx] = a2wide[1:1+ny,2:2+nx] + a2org
a2wide[2:2+ny,0:nx]   = a2wide[2:2+ny,0:nx] + a2org
a2wide[2:2+ny,1:1+nx] = a2wide[2:2+ny,1:1+nx] + a2org
a2wide[2:2+ny,2:2+nx] = a2wide[2:2+ny,2:2+nx] + a2org

a2wide[a2wide>0]=1
a2wide = a2wide[1:-1,1:-1]

#******* save ***********
datdir = "/home/utsumi/temp/bams2020/aphro"
sitepath = datdir + "/frac.site.%04d-%04d.npy"%(iY,eY)
widepath = datdir + "/mask.expand.site.th-%0.1f.%04d-%04d.npy"%(thratio, iY,eY)

np.save(sitepath, a2site.astype("float32"))
np.save(widepath, a2wide.astype("int16"))
print(sitepath)

###***** Figure *********
for idat in range(2):
    if figflag != True: continue

    if idat==0:
        a2fig = a2site
        #a2fig = ma.masked_less(a2site, thratio)
    else:
        a2fig = a2wide

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

    stitle = '%s %% %04d-%04d'%({0:"site",1:"expand"}[idat],iY,eY)

    axmap.set_title(stitle)

    #-- pcolormesh -------------------
    X,Y = np.meshgrid(a1lonbnd, a1latbnd)
    im = axmap.pcolormesh(X,Y, a2fig, cmap=cmap, norm=norm)
    plt.colorbar(im)
    plt.show()

    #plt.savefig(figpath)
    #print(figpath)


# %%