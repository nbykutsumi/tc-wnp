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
#figflag = False
figflag = True
difflag = True
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
#lens    = list(range(1,50+1))
#lens    = list(range(2,50+1))
lens    = list(range(1,1+1))
#lens    = list(range(5,50+1))

#lvname = ["steer"]
lvname = ["wind850"]
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

#************************************
# make grid distance data (for vorticity)
#************************************
"""
[x-1][x0][x1] --> calculate distance between x-1 to x1
"""
a1lat1 = a1latcnt
a1lat2 = a1latcnt
a1lon1 = np.full(len(a1lat1), a1loncnt[0])
a1lon2 = np.full(len(a1lat1), a1loncnt[1])
a1dist = util.calc_dist_gc_array(a1lat1,a1lat2,a1lon1,a1lon2)






#----------------------------------
for vname in lvname:
    #***************************************
    # Draw map (each Scen) 
    #***************************************
    for scen in lscen:
        if figflag != True: continue

        if vname=="steer":
            uvar = "usteer"
            vvar = "vsteer"

        elif vname=="wind850":
            uvar = "U850"
            vvar = "V850"

        a2varx= np.array([np.load( outbasedir + "/%s.%03d"%(scen,ens) + "/%s.%04d.%s.npy"%(uvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)
        a2vary = np.array([np.load( outbasedir + "/%s.%03d"%(scen,ens) + "/%s.%04d.%s.npy"%(vvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)
  
 
        X,Y = np.meshgrid(a1loncnt[x0:x1+1], a1latcnt[y0:y1+1])

        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
        axmap.coastlines()

        gl = axmap.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), xlocs=list(range(0,180,15)) )
        gl.ylabels_left  = True
        gl.ylabels_right = False
        gl.xlabels_top   = False
        gl.xlabels_bottom= True
    
        axmap.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    
        #-- coast lines ------
        axmap.coastlines()
    
        #-- vector ------------
        if vname=="steer":
            scale=1.5; width=0.2; minshaft=2
        elif vname=="wind850":
            scale=1.0; width=0.2; minshaft=2

        q = axmap.quiver(X, Y, a2varx, a2vary, angles="xy", units="xy", scale=scale, color="k", width=width, regrid_shape=15, minshaft=minshaft)


        #-- legend ------------
        axmap.quiverkey(q, 0.75, 0.95, 1, "1 m/s", labelpos="E", coordinates="figure")

        #-- title, figure name -----
        figpath = figdir + '/map.%s.%s.%04d-%04d.png'%(vname,scen,iY,eY)
        stitle = '%s %s\n'%(vname,scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])
    
        axmap.set_title(stitle)
    
        plt.show()
    
        plt.savefig(figpath)
        print(figpath)

    #***************************************
    # Draw map (difference) 
    #***************************************

    if vname=="steer":
        uvar = "usteer"
        vvar = "vsteer"

    a2hisx= np.array([np.load( outbasedir + "/%s.%03d"%("HPB",ens) + "/%s.%04d.%s.npy"%(uvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)
    a2hisy = np.array([np.load( outbasedir + "/%s.%03d"%("HPB",ens) + "/%s.%04d.%s.npy"%(vvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)

    a2natx= np.array([np.load( outbasedir + "/%s.%03d"%("HPB_NAT",ens) + "/%s.%04d.%s.npy"%(uvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)
    a2naty = np.array([np.load( outbasedir + "/%s.%03d"%("HPB_NAT",ens) + "/%s.%04d.%s.npy"%(vvar,Year,season))[y0:y1+1,x0:x1+1] for Year in lYear for ens in lens]).mean(axis=0)

    a2varx = a2hisx - a2natx
    a2vary = a2hisy - a2naty
 
    X,Y = np.meshgrid(a1loncnt[x0:x1+1], a1latcnt[y0:y1+1])

    fig = plt.figure(figsize=(6,4))
    axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
    axmap.coastlines()

    gl = axmap.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), xlocs=list(range(0,180,15)) )
    gl.ylabels_left  = True
    gl.ylabels_right = False
    gl.xlabels_top   = False
    gl.xlabels_bottom= True
    
    axmap.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    
    #-- coast lines ------
    axmap.coastlines()
    
    #-- speed differnce contour ------------
    if vname in ["steer"]:
        a2difspeed = np.sqrt(a2hisx**2+a2hisy**2) - np.sqrt(a2natx**2+a2naty**2) 

        #cmbnd = list(np.arange(-40,-5+0.01,5)) + [-2,2] + list(np.arange(5,40+0.01,5))
        cmbnd = list(np.arange(-1.1,1.1+0.01,0.2))
        cmlabels = list(map(float, ["%.1f"%x for x in cmbnd]))
        mycm = 'RdBu_r'

        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        Xbnd, Ybnd = np.meshgrid(a1lonbnd[x0:x1+2], a1latbnd[y0:y1+2])
        im  = plt.pcolormesh(Xbnd,Ybnd,a2difspeed, cmap=cmap, norm=norm)
        #-- draw colorbar ------
        cax = fig.add_axes([0.85,0.1,0.03,0.8])
        cbar= plt.colorbar(im, orientation='vertical', cax=cax)
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

    #-- vector ------------
    if vname=="steer":
        scale=0.3; width=0.2; minshaft=2
    elif vname=="wind850":
        scale=0.2; width=0.2; minshaft=2

    q = axmap.quiver(X, Y, a2varx, a2vary, angles="xy", units="xy", scale=scale, color="k", width=width, regrid_shape=15, minshaft=minshaft)

    #-- vector legend ------------
    axmap.quiverkey(q, 0.75, 0.95, 1, "1 m/s", labelpos="E", coordinates="figure")


    #-- title, figure name -----
    figpath = figdir + '/map.dif.%s.%04d-%04d.png'%(vname,iY,eY)
    stitle = 'dif %s\n'%(vname) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])
    
    axmap.set_title(stitle)
    
    plt.show()
    
    plt.savefig(figpath)
    print(figpath)



# %%
