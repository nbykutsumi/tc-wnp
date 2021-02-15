# %%
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
#%matplotlib inline
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

calcflag = True
#calcflag = False
figflag = True

#lieY = [[1951,1970]]
#lieY = [[1996,2015]]
lieY = [[1980,1997],[1998,2015]]

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


lrp = [1,5,10,20]
dpercent = {rp:100-100/(365*rp) for rp in lrp}
print(dpercent)


for (iY,eY) in lieY:
    if calcflag != True: continue

    lY  = range(iY,eY+1)
    lM  = range(1,12+1)
    
    ndays = (datetime(eY,12,31)-datetime(iY,1,1)).days +1
    
    a3dat = np.zeros([ndays, ny, nx], "float32")
    print(a3dat.shape)
    for Year in lY:
        print(Year)
        a3temp = aph.load_year(Year)[0]
        iday = (datetime(Year,1,1) - datetime(lY[0],1,1)).days
        a3dat[iday:iday+a3temp.shape[0],:,:] = a3temp
    
    
    for rp in lrp:
        percent = dpercent[rp]
        a2ptile = np.percentile(a3dat, percent, axis=0)
            
        #--- Save -----
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/aphro"
        util.mk_dir(exdir)
        expath= exdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iY, eY)
        np.save(expath, a2ptile)
        print(expath)
    
#***** Figure *********
for (iY,eY) in lieY:
    if figflag != True: continue

    for rp in lrp:
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/aphro"
        expath= exdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iY, eY)
        a2fig = np.load(expath) #  mm/day

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
        #cmbnd = list(np.arange(0,600,50))
        #cmlabels = list(np.arange(0,600,100))
        cmbnd = list(np.arange(0,300,20))
        cmlabels = list(np.arange(0,300,40))



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
        figdir  = "/home/utsumi/temp/bams2020/fig/map-prec"
        figpath = figdir + "/map.prec.rp-%03d.APHRO.%s.%04d-%04d.png"%(rp, regionfig, iY, eY)
        stitle = '%s-year R.P. mm/day APHRODITE\n'%(rp) + '%04d-%04d'%(iY,eY)

        axmap.set_title(stitle)

        #-- contourf -------------------
        X,Y = np.meshgrid(a1lonbnd, a1latbnd)
        #im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
        im = axmap.pcolormesh(X,Y, a2fig, cmap=cmap, norm=norm)
        plt.colorbar(im)
        plt.show()

        plt.savefig(figpath)
        print(figpath)




# %%
