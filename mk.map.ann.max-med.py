# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
#----------------------------------
import numpy as np
import scipy
import d4PDF
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from importlib import import_module
from numpy import ma

calcflag = False
figflag  = True
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lscen   = ['HPB']  # run={expr}-{scen}-{ens}
#lscen   = ['HPB_NAT'] # run={expr}-{scen}-{ens}
lscen   = ['HPB_NAT','HPB'] # run={expr}-{scen}-{ens}

lens    = list(range(1,50+1))
#lens    = list(range(21,50+1))
#lens    = [1,2]
miss_out = -9999.

iY,eY = 1990,2010
#iY,eY = 1990,1991
lY  = range(iY,eY+1)
ny,nx = 320, 640
hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'

detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)


for scen in lscen:
    maxbasedir = "/home/utsumi/temp/bams2020/max-prec-d"
    a3max = np.array([np.load(maxbasedir + "/%s.%03d/ann-max-prec.%04d.npy"%(scen,ens,Year)) for Year in lY for ens in lens])
    a2ave = a3max.mean(axis=0)
    a2med = np.median(a3max, axis=0)


    if calcflag ==True:
        maxdir = maxbasedir + "/%s.ave"%(scen)
        avepath= maxdir + "/ann-max-prec.ens-%03d-%03d.%04d-%04d.npy"%(lens[0],lens[-1],iY,eY)
        util.mk_dir(maxdir)  
        np.save(avepath, a2ave)
        print(avepath)

        maxdir = maxbasedir + "/%s.med"%(scen)
        medpath= maxdir + "/ann-max-prec.ens-%03d-%03d.%04d-%04d.npy"%(lens[0],lens[-1],iY,eY)
        util.mk_dir(maxdir)  
        np.save(medpath, a2med)
        print(medpath)

    print("-----------------------")
    if figflag != True: continue
    a2fig = a2med*60*60*24  # mm/day
    
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
    cmbnd = list(np.arange(0,200,10))
    cmlabels = list(np.arange(0,200,20))
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #-- title, figure name -----
    figdir  = '/home/utsumi/temp/bams2020/fig/map-prec'

    figpath = figdir + '/map.prec-max-d.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'median. annual max. precip. mm/day %s\n'%(scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 

    axmap.set_title(stitle)

    #-- contourf -------------------
    X,Y = np.meshgrid(a1lon, a1lat)
    im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
    plt.colorbar(im)
    plt.show()
                         
    plt.savefig(figpath)
    print(figpath)




# %%
