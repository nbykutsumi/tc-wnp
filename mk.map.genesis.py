# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
import socket
import glob
from bisect import bisect_left
#--------------------------------------
calcobj= True
#calcobj= False
#figobj = True
figobj = False

figmean = True
#figmean = False

ctype = 'TC'
#ctype = 'ALL'

iY, eY = 1990,2010
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])

#cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = np.arange(0,50,4)
#cmbnd = [0,1,5] + range(10,50,5)
#cmbnd = np.array(cmbnd)
#cmbnd = np.arange(0,50,4)
cmbnd = np.arange(0,50,4)
print(cmbnd)
#cmbnd = None
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen   = ["HPB","HPB_NAT"]
#lens    = range(1,2+1)
lens    = range(1,50+1)


detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
#IBTrACS     = import_module("%s.IBTrACS"%(detectName))


#wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outbaseDir = '/home/utsumi/temp/bams2020/map-genesis'
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-genesis'
util.mk_dir(figdir)
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
a1latinbnd = d4PDF.LatBnd()
a1loninbnd = d4PDF.LonBnd()
nyin    = len(a1latin)
nxin    = len(a1lonin)


dgrid = 5
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]  # for making count data, not for figure
#[[lllat,lllon],[urlat,urlon]] = [[10,100],[50,150]]  # for making count data, not for figure
a1latbnd = np.arange(lllat,urlat+0.01,dgrid)
a1lonbnd = np.arange(lllon,urlon+0.01,dgrid)
a1latcnt = (a1latbnd[:-1]+a1latbnd[1:])*0.5
a1loncnt = (a1lonbnd[:-1]+a1lonbnd[1:])*0.5

nysub = len(a1latcnt)
nxsub = len(a1loncnt)

miss_int= -9999

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

    if cmbnd is not None:
        extend = dpara['extend']
        cbar.set_ticks(cbar.ax.get_yticks())
        if extend=='both':
            cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:-1] + [""])

        if extend=='min':
            cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:])

        if extend=='max':
            cbar.set_ticklabels(list(cbar.ax.get_yticklabels())[:-1]+ [""])

        else:
            pass
    #-- coastline ---------------

    #-- Tiele -----
    stitle = dpara['title']
    axmap.set_title(stitle)
    #-- Save ------
    figpath = dpara['figpath']
    plt.savefig(figpath)
    plt.show()
    print(figpath)



#************************************
# d4PDF (Objective detection)
#************************************
thsstdeg = 27
exrvort = 3*1.0e-5
tcrvort = 3*1.0e-5
thwcore= 0
thdura = 36
thwind = 12
thwdif = -9999

thsstk     = thsstdeg +273.15
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsstdeg*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

for scen in lscen:
    for ens in lens:
        for Year in lYear:
            if calcobj is not True: continue

            a2num = np.zeros([nysub,nxsub],int16)

            for Mon in lMon:

                ibasedir = '/home/utsumi/temp/bams2020/tc-wise-WNP'
                idir = ibasedir + "/%s/%s.%03d/%04d%02d"%(slabel,scen,ens,Year,Mon)
                ssearch = idir + "/slp.*.bn"
                lslppath = sorted(glob.glob(ssearch))
                for slppath in lslppath:

                    #-- check duration ---
                    dura = os.path.getsize(slppath)  # float32
                    if dura < thdura: continue
                    #print(slppath)
                    srcdir= os.path.dirname(slppath)
                    fname = os.path.basename(slppath)
                    ipos = int(fname.split(".")[-2])

                    inity,initx = np.unravel_index(ipos, [nyin, nxin])

                    lat = a1latin[inity]
                    lon = a1lonin[initx]
                    ysub = int((lat - lllat)/dgrid)
                    xsub = int((lon - lllon)/dgrid)

                    if (ysub<0)or(ysub>=nysub)or(xsub<0)or(xsub>=nxsub): continue

                    a2num[ysub,xsub] = a2num[ysub,xsub] +1

            #-- Save --------
            print(a2num.sum())
            outdir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
            numpath= outdir + "/num.%04d.npy"%(Year)
            util.mk_dir(outdir)
            np.save(numpath, a2num.astype("int16"))
            print(numpath)
##***********************
## Figure: ensemble mean
##***********************
[[lllatfig,lllonfig],[urlatfig,urlonfig]] = [[10,105],[45,150]]


for scen in lscen:    
    outdir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)

    a3num = np.zeros([len(lens),nysub,nxsub],float32)

    for iens,ens in enumerate(lens):
        a2numtmp = np.array([np.load(outdir + "/num.%04d.npy"%(Year)).astype("int32") for Year in lYear]).sum(axis=0)
        a2numtmp = a2numtmp.astype("float32") / len(lYear) # count per year
    
        a3num[iens] = a2numtmp

    #--------------
    a2dat = a3num.mean(axis=0)
    print(a2dat.sum())
    dpara = {}
    dpara["cmbnd"] = np.arange(0,5+0.01,0.5)
    dpara["extend"]= "max"
    dpara["title"] = "TC genesis [#/year] %s"%(scen)
    dpara["figpath"]= figdir + "/map.genesis.%s.png"%(scen)
    draw_map(a2dat, dpara)
    


    #    dpara['title'] = 'Prob. of existence (d4PDF) [count/year] %04d-%04d'%(iY, eY) + '\n' + '%s-%s sst:%d ex:%.2f tc:%.2f \n wc:%.1f wind:%d wdif:%d ens-mean'%(expr, scen, thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif)
    #    dpara['figpath'] = figdir + '/map.freq.tc.obj.%s.%04d-%04d.ave.png'%(slabel, iY, eY)
    #    dpara['cmbnd'] = cmbnd
    #    dpara['extend']= "max"
    #    draw_map(a2fig, dpara)

# %%
