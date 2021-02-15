# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
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
calcflag = True
#calcflag = False
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
#lens    = list(range(1,3+1))
#lens    = list(range(1,50+1))
lens    = list(range(2,50+1))
#lens    = list(range(9,50+1))
#lens    = list(range(1,5+1))
#lens    = list(range(5,50+1))

#lvname = ["shear"]
#lvname = ["usteer","vsteer"]
#lvname = ["U850","V850"]
lvname = ["V850"]

region = "WNP"
#dBBox = {"WNP" :[[10,100],[50,250]]}
dBBox = {"WNP" :[[10,105],[45,150]]}
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

d4 = d4PDF.snp_6hr_2byte(vtype="atm",dbbaseDir=d4pdfdir)
a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)

ny = d4.ny
nx = d4.nx
print(type(ny))

[[lllat,lllon],[urlat,urlon]] = dBBox[region]
a1latcnt = d4PDF.Lat()
a1loncnt = d4PDF.Lon()

a1latbnd = d4PDF.LatBnd()
a1lonbnd = d4PDF.LonBnd()
y0 = bisect_left(a1latbnd, lllat)
y1 = bisect_left(a1latbnd, urlat)
x0 = bisect_left(a1lonbnd, lllon)
x1 = bisect_left(a1lonbnd, urlon)

#----------------------------------
#************************************
# d4PDF (Objective detection)
#************************************
for vname in lvname:
    for scen in lscen:
        if calcflag !=True: continue
    
        for ens in lens:
            for Year in lYear:
                a3var = np.zeros([len(lMon),ny,nx],"float32")
                for i,Mon in enumerate(lMon):
                    print(vname,scen,ens,Year,Mon)
                    if vname in ["usteer","vsteer"]:
                        vname300 = "%s300"%(str.upper(vname[0]))
                        vname500 = "%s500"%(str.upper(vname[0]))
                        vname850 = "%s850"%(str.upper(vname[0]))

                        alw = 850/float(300+500+850)
                        amd = 500/float(300+500+850)
                        aup = 300/float(300+500+850)

                        a2lw = d4.load_ave_mon(vname850, scen, ens, Year, Mon)
                        a2md = d4.load_ave_mon(vname500, scen, ens, Year, Mon)
                        a2up = d4.load_ave_mon(vname300, scen, ens, Year, Mon)
                        print(a2up.shape,a3var.shape)
                        a2vartmp= a2lw*alw + a2md*amd + a2up*aup

                    if vname == "shear":
                        a2ulw = d4.load_ave_mon("U850", scen, ens, Year, Mon)
                        a2uup = d4.load_ave_mon("U300", scen, ens, Year, Mon)

                        a2vlw = d4.load_ave_mon("V850", scen, ens, Year, Mon)
                        a2vup = d4.load_ave_mon("V300", scen, ens, Year, Mon)

                        a2du = a2uup - a2ulw
                        a2dv = a2vup - a2vlw
                        a2vartmp = np.sqrt(a2du**2 + a2dv**2)
                    else:
                        a2vartmp = d4.load_ave_mon(vname, scen, ens, Year, Mon)

                    a3var[i] = a2vartmp.astype("float32")

                a2var = ma.masked_invalid(a3var).mean(axis=0).filled(np.nan)
                #-- save ---
                outdir = outbasedir + "/%s.%03d"%(scen,ens)
                varpath = outdir + "/%s.%04d.%s.npy"%(vname,Year,season)
                util.mk_dir(outdir)
                np.save(varpath, a2var.astype("float32"))
                print(varpath)

    ##*** Draw map (each Scen) *************
    for scen in lscen:
        if figflag !=True: continue

        a3var = np.array([np.load( outbasedir + "/%s.%03d"%(scen,ens) + "/%s.%04d.%s.npy"%(vname,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens])
   
        a2fig = a3var.mean(axis=0)

        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
        axmap.coastlines()

        gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        xticks = np.arange(-180,180+1,15)
        yticks = np.arange(-90,90+1,15)
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
        axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

        axmap.set_extent([lllon,urlon,lllat,urlat])

        #-- Make new colormap (white at the lower end) --
        upper = matplotlib.cm.jet(np.arange(256))
        #lower = np.ones((int(256/4),4))
        lower = np.ones((int(256/8),4))
        for i in range(3):
            lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
        mycm = np.vstack(( lower, upper ))
        mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

        #-- color boundaries norm ------
        if vname=="shear":
            cmbnd = list(np.arange(0,30,2))
            cmlabels = list(np.arange(0,30,4))


        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        #--coast lines ------
        axmap.coastlines()

        #--pcolormesh  ----------
        Xbnd, Ybnd = np.meshgrid(a1lonbnd[x0:x1+2], a1latbnd[y0:y1+2])
        im = axmap.pcolormesh(Xbnd,Ybnd, a2fig, cmap=mycm, norm=norm)

        #-- colorbar -------------
        cbar = plt.colorbar(im)
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

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
        if vname in ["usteer","vsteer","U850","V850"]: continue
        a3his = np.array([np.load( outbasedir + "/%s.%03d"%("HPB",ens) + "/%s.%04d.%s.npy"%(vname,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens])

        a3nat = np.array([np.load( outbasedir + "/%s.%03d"%("HPB_NAT",ens) + "/%s.%04d.%s.npy"%(vname,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens])

        a2tv, a2pv = scipy.stats.ttest_ind(a3his, a3nat, axis=0, equal_var=False, nan_policy="omit")
        a2fig = a3his.mean(axis=0) - a3nat.mean(axis=0)

        a2hatch = ma.masked_where(a2pv>0.05, a2fig)

        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
        axmap.coastlines()

        gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        xticks = np.arange(-180,180+1,15)
        yticks = np.arange(-90,90+1,15)
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
        axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

        axmap.set_extent([lllon,urlon,lllat,urlat])

        #-- color boundaries norm --------
        if vname=="shear":
            cmbnd = list(np.arange(-1.7,1.7+0.01,0.2))
            cmlabels = list(map(float, ["%.1f"%x for x in cmbnd]))

        mycm = 'RdBu_r'
        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        #-- coast lines ------
        axmap.coastlines()
    
        #--pcolormesh ----------
        Xbnd, Ybnd = np.meshgrid(a1lonbnd[x0:x1+2], a1latbnd[y0:y1+2])
        im = axmap.pcolormesh(Xbnd,Ybnd, a2fig, cmap=cmap, norm=norm)

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
