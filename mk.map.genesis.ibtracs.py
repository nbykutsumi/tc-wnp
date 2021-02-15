# %%
import matplotlib
%matplotlib inline
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import IBTrACS
import myfunc.util as util
from datetime import datetime, timedelta
import os, sys
from bisect import bisect_left
#-------------------------------
calcflag = False

iYear = 1980
eYear = 2010
lYear = list(range(iYear,eYear+1))
ib = IBTrACS.IBTrACS()

a1latbnd = np.arange(-90,90+0.1,5)
a1lonbnd = np.arange(0,360+0.1,5)
lat0  = -90
lon0  = 0
dgrid = 5
ny = len(a1latbnd)-1
nx = len(a1lonbnd)-1

#**************************
def draw_map(a2dat, dpara):
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())

    gl        = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,901, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)

    axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
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

    cmlabels = dpara['cmlabels']
    if cmlabels is not None:
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

    #-- coastline ---------------

    #-- Tiele -----
    stitle = dpara['title']
    axmap.set_title(stitle)
    #-- Save ------
    figpath = dpara['figpath']
    plt.savefig(figpath)
    plt.show()
    print(figpath)


#********************************************************
for Year in lYear:
    if calcflag !=True: continue

    idtime = datetime(Year,1,1,0)
    edtime = datetime(Year,12,31,18)
    dtc = ib.ret_dlonlat(idtime,edtime, dhours=6, lvaridx=[0])
    ldtime= list(dtc.keys())
    ltc = list(dtc.values())
    
    ltc = list(zip(ldtime, ltc))
    ltc = sorted(ltc, key=lambda x:x[0]) # sort by datetime
    ltc = list(zip(*ltc))[1]  # take only values
    ltc = [x for l in ltc for x in l]
    
    a1lon, a1lat, a1tcid = list(map(np.array, list(zip(*ltc))))
    atcid = []
    alon  = []
    alat  = []
    for (lon,lat,tcid) in zip(a1lon, a1lat, a1tcid):
        if tcid not in atcid: 
            atcid.append(tcid)
            alon.append(lon) 
            alat.append(lat) 
 
    atcid = np.array(atcid) 
    alon  = np.array(alon)
    alat  = np.array(alat)

    a1x = ((alon - lon0)/dgrid).astype(int16)
    a1y = ((alat - lat0)/dgrid).astype(int16)

    a2num = np.zeros([ny,nx], "int16")
    for i in range(len(a1x)):
        x=a1x[i]
        y=a1y[i]
        a2num[y,x] = a2num[y,x]+1

    #-- save ----
    datdir = '/home/utsumi/temp/bams2020/map-genesis/ibtracs'
    util.mk_dir(datdir)
    datpath= datdir + "/num.%04d.npy"%(Year)
    np.save(datpath, a2num)
    print(datpath) 


#******************************************
# Figure
#******************************************
[[lllatfig,lllonfig],[urlatfig,urlonfig]] = [[10,105],[45,150]]

datdir = '/home/utsumi/temp/bams2020/map-genesis/ibtracs'
a3num = np.array([np.load(datdir + "/num.%04d.npy"%(Year)) for Year in lYear])

a2dat = a3num.mean(axis=0)
print(a2dat.sum())

figdir  = '/home/utsumi/temp/bams2020/fig/map-genesis'
util.mk_dir(figdir)

dpara = {}
#dpara["cmbnd"] = np.arange(0,5+0.01,0.5)
dpara["cmbnd"] = [0,0.1,0.5,1,1.5,2,2.5,3,3.5,4]
dpara["cmlabels"] = [0,0.1,0.5,1,1.5,2,2.5,3,3.5,4]
dpara["extend"]= ""
dpara["title"] = "TC genesis [#/year] Best-track"
dpara["figpath"]= figdir + "/map.genesis.ibtracs.png"
draw_map(a2dat, dpara)

x0 = bisect_left(a1lonbnd, lllonfig)
x1 = bisect_left(a1lonbnd, urlonfig)
y0 = bisect_left(a1latbnd, lllatfig)
y1 = bisect_left(a1latbnd, urlatfig)

a2dat[y0:y1+1,x0:x1+1].sum()
# %%
