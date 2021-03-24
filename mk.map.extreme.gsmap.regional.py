# %%
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker

import numpy as np
import scipy
import d4PDF
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from bisect import bisect_left
from importlib import import_module
from numpy import ma
import GSMaP
import d4PDF
import myfunc.regrid.Regrid as Regrid
#calcflag = True
calcflag = False
figflag = True

#regrid = True
regrid = False

if regrid==True:
    sregrid="UP"
else:
    sregrid="ORG"
    
iY,eY = 2001,2010
#iY,eY = 2001,2001
#iY,eY = 2009,2010

iDTime = datetime(iY,1,1)
eDTime = datetime(eY,12,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

ny,nx = 320, 640
hostname = socket.gethostname()
#if hostname =='shui':
#    dbDir    = "/media/disk2/data/GSMaP" 
#elif hostname=='well':
#    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

gs = GSMaP.GSMaP_daily(prj="standard",ver="v6",prdName="daily_G", compressed=True)

a1latorg = gs.Lat   # 0.1 
a1lonorg = gs.Lon   # 0.1

if regrid==True:
    a1latup = d4PDF.Lat()   # dlat ~ 0.5615674
    a1lonup = d4PDF.Lon()   # dlon = 0.5625
else:
    a1latup = a1latorg
    a1lonup = a1lonorg


region= "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

y0 = bisect_left(a1latup, lllat)
y1 = bisect_left(a1latup, urlat)
x0 = bisect_left(a1lonup, lllon)
x1 = bisect_left(a1lonup, urlon)

nyreg = y1-y0+1    # 
nxreg = x1-x0+1


lrp = [1,5,10,20]
dpercent = {rp:100-100/(365*rp) for rp in lrp}


us = Regrid.UpScale()
us(a1latorg, a1lonorg, a1latup, a1lonup, globflag=True)

ndays = len(lDTime)
if calcflag == True:
    a3dat = np.zeros([ndays, nyreg, nxreg], "float32")
    for iday,dtime in enumerate(lDTime):
        print(dtime)
        a2org = gs.load_day(dtime)

        if regrid==True:
            a2up  = us.upscale(a2org, pergrid=False, miss_in=-9.9990002e+2, miss_out=-9999)
        else:
            a2up  = a2org

        a3dat[iday,:,:] = a2up[y0:y1+1,x0:x1+1]
 
    for rp in lrp:
        percent = dpercent[rp]
        a2ptile = np.percentile(a3dat, percent, axis=0)
    
        print(rp)
        print(a2ptile.max())
        #--- Save -----
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/GSMaP.%s.%04d-%04d"%(sregrid,iY,eY)
        util.mk_dir(exdir)
        expath= exdir + "/prec.rp-%03d.%s.npy"%(rp, region)
        np.save(expath, a2ptile)
        print(expath)


#***** Figure *********
for rp in lrp:
    if figflag != True: continue
    exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/GSMaP.%s.%04d-%04d"%(sregrid,iY,eY)
    expath= exdir + "/prec.rp-%03d.%s.npy"%(rp, region)
    a2fig = np.load(expath) * 24  # mm/h --> mm/day
    print(expath)
    print(a2fig.max())

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
    #clevs = range(-20,20+1,4)
    cmbnd = [0,40,80,120,160,200,300,400,500]

    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    dbboxfig = {"WNP":[[0,100],[50,150]]}
    [[lllatfig,lllonfig],[urlatfig,urlonfig]] = dbboxfig[region]
    axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
    axmap.coastlines()

    #-- title, figure name -----
    figdir  = "/home/utsumi/temp/bams2020/fig/map-prec"
    figpath = figdir + "/map.prec.rp-%03d.GSMaP.%s.%s.%04d-%04d.png"%(rp, sregrid, region,iY,eY)
    stitle = '%s-year R.P. mm/day GSMaP %s\n'%(rp, sregrid) + '%04d-%04d'%(iY,eY)

    axmap.set_title(stitle)

    #-- contourf -------------------
    a1latreg = a1latup[y0:y1+1]
    a1lonreg = a1lonup[x0:x1+1]

    X,Y = np.meshgrid(a1lonreg, a1latreg)
    im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
    plt.colorbar(im)
    plt.show()

    plt.savefig(figpath)
    print(figpath)


# %%
