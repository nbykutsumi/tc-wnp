# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
from bisect import bisect_left
from detect_fsub import *
import socket
#--------------------------------------
#calcflag = True
calcflag = False
figflag  = True
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
lens    = [1,5]
noleap  = False

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
figdir  = '/home/utsumi/temp/bams2020'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640
miss =-9999
miss_int= -9999

#************************************
# d4PDF 
#************************************
miss =-9999
for scen in lscen:
    if calcflag != True: continue

    for ens in lens:
        print(('ens=',ens))
        for Year in lYear:
            ldtime = util.ret_lDTime(datetime(Year,1,1,0), datetime(Year,12,31,18), timedelta(hours=6))

            a2prec = np.array([d4sfc.load_6hr_mon("PRECIPI", scen, ens, Year, Mon).sum(axis=0) for Mon in lMon]).sum(axis=0) / float(len(ldtime))

            #--- Save file ---------
            precbasedir = "/home/utsumi/temp/bams2020/mqp-prec-all"
            precdir = precbasedir + "/%s.%03d"%(scen,ens)
            util.mk_dir(precdir)

            precpath = precdir + "/prec.%04d.npy"%(Year)
            np.save(precpath, a2prec.astype("float64"))
            print(precpath)

###*** Draw map *************
for scen in lscen:
    if figflag !=True: continue

    a3prec = np.empty([len(lens), ny,nx], "float64")
    for iens, ens in enumerate(lens):
        precbasedir = "/home/utsumi/temp/bams2020/mqp-prec-all"
        a3prec[iens] = np.array([np.load(precbasedir + "/%s.%03d/prec.%04d.npy"%(scen,ens,Year)) for Year in lYear]).mean(axis=0)

    #if scen == "HPB":
    #    a3prec_his = a3prec

    #elif scen=="HPB_NAT":
    #    a3prec_nat = a3prec

    a2fig = a3prec.mean(axis=0)*60*60*24 *365
    
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
    #clevs = range(-20,20+1,4)
    cmbnd = list(np.arange(0,4000,200))
    cmlabels = list(np.arange(-11,11+1,4))
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #-- Tiele -----

    #-- title, figure name -----
    figpath = figdir + '/map.prec-all.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'precipitation mm/year %s\n'%(scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 

    axmap.set_title(stitle)

    #-- contourf -------------------
    X,Y = np.meshgrid(a1lon, a1lat)
    im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
    plt.colorbar(im)
    plt.show()
                         
    plt.savefig(figpath)
    print(figpath)
