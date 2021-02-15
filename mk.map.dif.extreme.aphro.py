# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import scipy.stats
#----------------------------------
import sys, os, pickle
#from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import util
import calendar
import socket
import APHRODITE
#--------------------------------------
lys = [
    [[1951,1970],[1996,2015]]
]
#-----------------
figdir  = '/home/utsumi/temp/bams2020/fig/map-prec'
util.mk_dir(figdir)

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


lrp = [1,5,10,20]
dpercent = {rp:100-100/(365*rp) for rp in lrp}
print(dpercent)

lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
[[lllat,lllon],[urlat,urlon]] = [[10,105],[45,150]] # for figure

#************************************
for ([iY0,eY0],[iY1,eY1]) in lys:
    for rp in lrp:
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/aphro"

        expath0= exdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iY0, eY0)
        expath1= exdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iY1, eY1)

        a2dat0 = np.load(expath0)
        a2dat1 = np.load(expath1)
        a2dif = a2dat1 - a2dat0   # mm/day

        #***********************
        # Figure:
        #***********************
        a2fig = a2dif

        #-- title, figure name -----
        stitle = 'diff. %d-year R.P mm/day APHRODITE\n(%04d-%04d to %04d-%04d)'%(rp, iY0,eY0,iY1,eY1) 
        figpath = figdir + '/map.dif.prec.rp-%03d.APHRO.%04d-%04d.%04d-%04d.png'%(rp,iY0,eY0,iY1,eY1)
        #---------------------------
        figmap   = plt.figure(figsize=(6,4))
        axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())

        #-- grid lines ---
        gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        xticks   = np.arange(-180, 180+1, 15)
        yticks   = np.arange(-90,90, 15)
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        axmap.set_xticks(xticks, crs=ccrs.PlateCarree())
        axmap.set_yticks(yticks, crs=ccrs.PlateCarree())
        #-- set extent and coastlines----
        axmap.set_extent([lllon,urlon,lllat,urlat])
        axmap.coastlines(color="k")

        #-- color boundaries norm --------
        cmbnd = list(range(-90,90+1,20))
        cmlabels = list(range(-90,90+1,20))
        mycm = 'RdBu_r'

        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        #-- contourf --------------
        X,Y = np.meshgrid(a1lonbnd,a1latbnd)
        #vmin, vmax = cmbnd[0], cmbnd[-1]
        #im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
        im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm)
        print(cmbnd)

        #-- draw colorbar ------
        cax = figmap.add_axes([0.84,0.2,0.03,0.6])
        cbar= plt.colorbar(im, orientation='vertical', cax=cax)
        #cbar.set_ticks(cbar.ax.get_yticks())
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

        #-- Tiele -----
        axmap.set_title(stitle)

        #-- Save ------
        plt.savefig(figpath)
        plt.show()
        print(figpath)


        # %%
