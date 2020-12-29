# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#----------------------------------
import sys, os, pickle
from   numpy import *
import scipy.stats
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
from bisect import bisect_left
from detect_fsub import *
import socket
#--------------------------------------
#iY = 1990
#eY = 2010

iY = 1990
eY = 2010

lYear = range(iY,eY+1)
lMon  = range(1,12+1)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
#lens    = list(range(16,20+1))
#lens    = list(range(22,50+1))
#lens    = list(range(36,50+1))
#lens    = list(range(1,50+1))
#lens    = list(range(21,50+1))
lens    = list(range(1,5+1))
region  = "WNP"
detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'


compbasedir= '/home/utsumi/temp/bams2020/composite'
figdir  = '/home/utsumi/temp/bams2020/fig/map-prec'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640
miss =-9999
miss_int= -9999

#************************************
for scen in lscen:
    precbasedir  = "/home/utsumi/temp/bams2020/tc-prec-all"
    bboxpath = precbasedir + "/HPB.001/bbox.%s.npy"%(region)
    yxbboxpath = precbasedir + "/HPB.001/yxbbox.%s.npy"%(region)
    bbox   = np.load(bboxpath)
    yxbbox = np.load(yxbboxpath)

    [[lllat,lllon],[urlat,urlon]] = bbox
    [[y0,x0],[y1,x1]] = yxbbox
    nyreg = y1-y0+1
    nxreg = x1-x0+1

    a3max = np.empty([len(lens), nyreg,nxreg], "float64")
    for iens, ens in enumerate(lens):

        a3tmp = np.empty([len(lYear),nyreg,nxreg], "float64")
        for itmp,Year in enumerate(lYear):
            precpath = precbasedir + "/%s.%03d/prec-tc.%s.%03d.npy"%(scen,ens,region,Year)
        
            a3prec = np.load(precpath)

            a3tmp[itmp]= a3prec.max(axis=0)

        a3max[iens] = a3tmp.max(axis=0)

    #if scen == "HPB":
    #    a3dat_his = a3max

    #elif scen=="HPB_NAT":
    #    a3dat_nat = a3max


    #***** Figure *********
    a2fig = np.median(a3max, axis=0)*60*60*24  # mm/day

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
    cmbnd = list(np.arange(0,600,50))
    cmlabels = list(np.arange(0,600,100))
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    dbboxfig = {"WNP":[[0,100],[50,150]]}
    [[lllatfig,lllonfig],[urlatfig,urlonfig]] = dbboxfig[region]
    axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
    axmap.coastlines()

    #-- Tiele -----

    #-- title, figure name -----
    figpath = figdir + '/map.prec-tc-max.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'median of max tc prec. mm/year %s\n'%(scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 

    axmap.set_title(stitle)

    #-- contourf -------------------
    a1latreg = a1lat[y0:y1+1]
    a1lonreg = a1lon[x0:x1+1]

    X,Y = np.meshgrid(a1lonreg, a1latreg)
    im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
    plt.colorbar(im)
    plt.show()
                         
    plt.savefig(figpath)
    print(figpath)




# %%
