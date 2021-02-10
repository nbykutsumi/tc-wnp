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
import scipy.stats
#--------------------------------------
#calcflag = True
calcflag = False
#figflag  = True
figflag  = False
difflag  = True
iY = 1990
eY = 2010
lYear = range(iY,eY+1)

season = "MJJASO"
lMon  = {"ALL":range(1,12+1),
         "MJJASO":range(5,10+1)}[season]

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
lens    = list(range(1,50+1))
#lens    = list(range(2,50+1))
#lens    = list(range(1,2))

#lvname = ["TGEF","SLP"]
lvname = ["TGEF"]
#lvname = ["SLP"]

region = "NP"
dBBox = {"NP" :[[0,100],[50,250]],
         "WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'


outbasedir = '/home/utsumi/temp/bams2020/map-env'
figdir  = '/home/utsumi/temp/bams2020/fig/map-env'
util.mk_dir(figdir)

d4 = d4PDF.avr_mon_320x640(vtype="sfc", dbbaseDir=d4pdfdir)
a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
#************************************
# d4PDF (Objective detection)
#************************************
for vname in lvname:
    for scen in lscen:
        if calcflag !=True: continue
    
        for ens in lens:
            print(vname,scen,ens)
            a2var = np.array([d4.load_ave_mon(vname, scen, ens, Year, Mon) for Mon in lMon for Year in lYear]).mean(axis=0)

            #-- save ---
            outdir = outbasedir + "/%s.%03d"%(scen,ens)
            varpath = outdir + "/%s.%04d-%04d.%s.npy"%(vname,iY,eY,season)
            util.mk_dir(outdir)
            np.save(varpath, a2var.astype("float32"))
            print(varpath)

    ##*** Draw map (each Scen) *************
    for scen in lscen:
        if figflag !=True: continue

        [[lllat,lllon],[urlat,urlon]] = dBBox[region]
        a1lat = d4PDF.Lat()
        a1lon = d4PDF.Lon()

        a3var = np.array([np.load( outbasedir + "/%s.%03d"%(scen,ens) + "/%s.%04d-%04d.%s.npy"%(vname,iY,eY,season)) for ens in lens])
   
        a2fig = a3var.mean(axis=0)
        if   vname == "TGEF":
            a2fig = a2fig -273.15 # K --> deg.C
        elif vname == "SLP": 
            a2fig = a2fig * 0.01  # Pa --> hPa
    
        X,Y = np.meshgrid(a1lon, a1lat)

        clon= 180
        fig = plt.figure(figsize=(6,2))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree(central_longitude=clon))
        axmap.coastlines()

        gl = axmap.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), xlocs=[100,120,140,160,180,-160,-140,-120,-100]) 
        gl.ylabels_left  = True
        gl.ylabels_right = False
        gl.xlabels_top   = False
        gl.xlabels_bottom= True
    
        axmap.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    
        #-- Make new colormap (white at the lower end) --
        upper = matplotlib.cm.jet(np.arange(256))
        #lower = np.ones((int(256/4),4))
        lower = np.ones((int(256/8),4))
        for i in range(3):
            lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
        mycm = np.vstack(( lower, upper ))
        mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])
    
        #-- color boundaries norm ------
        if vname=="TGEF":
            cmbnd = list(np.arange(5,33+1,2))
            #cmlabels = list(np.arange(0,50+1,10))
        elif vname=="SLP":
            cmbnd = list(np.arange(1004,1024+1,1))

        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        #-- coast lines ------
        axmap.coastlines()
    
        #--contourf ----------
        im = axmap.contourf(X-clon,Y, a2fig, levels=cmbnd, cmap=mycm, extend="max")
        plt.colorbar(im)
    
        #-- title, figure name -----
        figpath = figdir + '/map.%s.%s.%04d-%04d.png'%(vname,scen,iY,eY)
        stitle = '%s %s\n'%(vname,scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])
    
        axmap.set_title(stitle)
    
        plt.show()
    
        plt.savefig(figpath)
        print(figpath)

    #************************************
    # Draw map (difference) 
    #************************************
    if difflag ==True:
        [[lllat,lllon],[urlat,urlon]] = dBBox[region]
        a1lat = d4PDF.Lat()
        a1lon = d4PDF.Lon()
        a1latbnd = d4PDF.LatBnd()
        a1lonbnd = d4PDF.LonBnd()

        a3his = np.array([np.load( outbasedir + "/%s.%03d"%("HPB",ens) + "/%s.%04d-%04d.%s.npy"%(vname,iY,eY,season)) for ens in lens])

        a3nat = np.array([np.load( outbasedir + "/%s.%03d"%("HPB_NAT",ens) + "/%s.%04d-%04d.%s.npy"%(vname,iY,eY,season)) for ens in lens])


        a2tv, a2pv = scipy.stats.ttest_ind(a3his, a3nat, axis=0, equal_var=False, nan_policy="omit")
        a2fig = a3his.mean(axis=0) - a3nat.mean(axis=0)

        if vname == "SLP": 
            a2fig = a2fig * 0.01  # Pa --> hPa
    
        a2hatch = ma.masked_where(a2pv>0.05, a2fig)
        X,Y = np.meshgrid(a1lon, a1lat)

        clon= 180
        fig = plt.figure(figsize=(6,2))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree(central_longitude=clon))
        axmap.coastlines()

        gl = axmap.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), xlocs=[100,120,140,160,180,-160,-140,-120,-100]) 
        gl.ylabels_left  = True
        gl.ylabels_right = False
        gl.xlabels_top   = False
        gl.xlabels_bottom= True
    
        axmap.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    
        #-- color boundaries norm --------
        if vname=="TGEF":
            #cmbnd = list(np.arange(-1.9,1.9+0.01,0.2))
            cmbnd = list(map(float, ["%.1f"%x for x in cmbnd]))
            cmlabels = [-1.9,-1.5,-1.1,-0.7,-0.3, 0.3, 0.7, 1.1, 1.5, 1.9]
        elif vname=="SLP":
            cmbnd = list(np.arange(-0.9,0.9+0.01,0.2))
            cmlabels = cmbnd

        mycm = 'RdBu_r'
        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        #-- coast lines ------
        axmap.coastlines()
    
        #--contourf ----------
        im = axmap.contourf(X-clon,Y, a2fig, levels=cmbnd, cmap=mycm, extend="both")
   
        #-- colorbar ---
        cbar = plt.colorbar(im)
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

        #-- hatch ------------
        #Xbnd,Ybnd = a1lonbnd, a1latbnd
        #plt.pcolor(Xbnd-clon, Ybnd, a2hatch, hatch="/", alpha=0.)

        #-- title, figure name -----
        figpath = figdir + '/map.dif.%s.%04d-%04d.png'%(vname,iY,eY)
        stitle = 'diff. %s\n'%(vname) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])
    
        axmap.set_title(stitle)
    
        plt.show()
    
        plt.savefig(figpath)
        print(figpath)









# %%
