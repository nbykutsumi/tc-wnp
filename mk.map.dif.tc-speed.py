# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

#----------------------------------
import sys, os, pickle, socket
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
from bisect import bisect_left
from collections import deque
from scipy import stats
import scipy
from math import sin, cos, acos

#--------------------------------------
figflag = True

#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 1990,1991
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}

dlatsub  = 5
dlonsub  = 5
llatsub0 = np.arange(5,40+0.1,dlatsub)[::-1]
llonsub0 = np.arange(105,145+0.1,dlonsub)
lsubbox = [[[latsub0,lonsub0],[latsub0+dlatsub,lonsub0+dlonsub]] for latsub0 in llatsub0 for lonsub0 in llonsub0]

[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]   # for loading tc-vector data

detectName = 'wsd_d4pdf_20201209-py38'
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))  # test
d4PDF       = import_module("%s.d4PDF"%(detectName))
#IBTrACS     = import_module("%s.IBTrACS"%(detectName))

hostname=socket.gethostname()
if hostname=="shui":
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
elif hostname=="well":
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
else:
    print("check hostname",hostname)
    sys.exit()

vectbaseDir = '/home/utsumi/temp/bams2020/vect-tc-var'
speedbasedir= "/home/utsumi/temp/bams2020/map-speed"
figdir  = '/home/utsumi/temp/bams2020/fig/speed'
util.mk_dir(figdir)
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

dgrid = 5.0
#dgrid = 2.5
o1latbnd= np.arange(-90,90+0.01,dgrid)
o1lonbnd= np.arange(0,360+0.01,dgrid)
o1latcnt= np.arange(-90+dgrid*0.5, 90+0.01,dgrid)
o1loncnt= np.arange(0+dgrid*0.5, 360+0.01,dgrid)

nyout = len(o1latcnt)
nxout = len(o1loncnt)

miss_int= -9999

#************************************
# d4PDF (Objective detection)
#************************************
thsst   = 27
exrvort = 3*1.0e-5
tcrvort = 3*1.0e-5
thwcore= 0
thdura = 36
thwind = 12
thwdif = -9999
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5
slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

#-- Draw ----
for scen in lscen:
    speeddir = speedbasedir + "/%s"%(slabel)
    avepath= speeddir + "/ave.%s.%04d-%04d.npy"%(scen,iY,eY)  # contains 50 ensemble
    sumpath= speeddir + "/sum.%s.%04d-%04d.npy"%(scen,iY,eY)

    sumxpath= speeddir + "/sumx.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumypath= speeddir + "/sumy.%s.%04d-%04d.npy"%(scen,iY,eY)
    numpath= speeddir + "/num.%s.%04d-%04d.npy"%(scen,iY,eY)

    a3sum  = ma.masked_invalid(np.load(sumpath))
    a3sumx = ma.masked_invalid(np.load(sumxpath))
    a3sumy = ma.masked_invalid(np.load(sumypath))
    a3num  = ma.masked_invalid(np.load(numpath))

    #a3ave  = np.load(avepath)
    a3ave  = ma.masked_where(a3num==0, a3sum) / a3num
    a3ave  = ma.masked_invalid(a3ave).filled(np.nan)

    a2sum = a3sum.sum(axis=0)
    a2num = a3num.sum(axis=0)
    a2ave = ma.masked_where(a2num==0, a2sum) / a2num

    a2sumx= a3sumx.sum(axis=0)
    a2sumy= a3sumy.sum(axis=0)
    a2avex= ma.masked_where(a2num==0, a2sumx) / a2num
    a2avey= ma.masked_where(a2num==0, a2sumy) / a2num

    if scen=="HPB":
        a3avehis = a3ave
        a2avehis = a2ave
        a2avexhis = a2avex
        a2aveyhis = a2avey
    elif scen=="HPB_NAT":
        a3avenat = a3ave
        a2avenat = a2ave
        a2avexnat = a2avex
        a2aveynat = a2avey
    else:
        print("check scen",scen)
        sys.exit()

#***********************
#***********************
a2dif = a2avehis - a2avenat
a2difx = a2avexhis - a2avexnat
a2dify = a2aveyhis - a2aveynat

a2difrat = (ma.masked_where(a2avenat==0, a2avehis - a2avenat)/a2avenat).filled(np.nan) * 100


#if ma.is_masked(a3avehis): a3avehis=a3avehis.filled(np.nan)
#if ma.is_masked(a3avenat): a3avenat=a3avenat.filled(np.nan)

a2tv, a2pv = scipy.stats.ttest_ind(a3avehis, a3avenat, axis=0, equal_var=False, nan_policy="omit")

#****************************************************
# Figure: difference in km/h
#****************************************************
a2fig = a2dif
a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#a2hatch = ma.masked_where(a2fig<0., a2fig)

#---------------------------
figmap   = plt.figure(figsize=(6,4))
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())

#-- grid lines ---
gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,90, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

#-- set extent and coastlines----
#axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.set_extent([105,150,10,45])  # for figure
axmap.coastlines(color="k")

#-- color boundaries norm --------
cmbnd = list(np.arange(-5,5+0.1,1))
cmlabels = list(np.arange(-5,5+0.1,1))

#cmbnd = list(np.arange(0,1+0.001,0.05))   # test
#cmlabels = list(np.arange(0,1+0.001,0.05))  # test

mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

#-- contourf --------------
X,Y = np.meshgrid(o1loncnt,o1latcnt)
#vmin, vmax = cmbnd[0], cmbnd[-1]
im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
print(cmbnd)

#-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -1, 1) # should be adjusted considering
plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)

#-- draw colorbar ------
cax = figmap.add_axes([0.85,0.1,0.03,0.7])
cbar= plt.colorbar(im, orientation='vertical', cax=cax)
cbar.set_ticks(cmlabels)
cbar.set_ticklabels(cmlabels)

#-- vector ------------
#axmap.quiver(X, Y, a2difx, a2dify, angles="xy", units="xy", scale=1, color="gray", width=0.4)

#-- title, figure name -----
nens = a3avehis.shape[0]
figpath = figdir + '/map.dif.speed.%04d-%04d.png'%(iY,eY)
stitle = 'TC speed diff. (km/h)\n' + '%04d-%04d nens:%03d'%(iY,eY,nens)

axmap.set_title(stitle)

plt.show()

plt.savefig(figpath)
print(figpath)


#****************************************************
# Figure: difference in ratio (% increase)
#****************************************************
a2fig = a2difrat
a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#a2hatch = ma.masked_where(a2fig<0., a2fig)

#---------------------------
figmap   = plt.figure(figsize=(6,4))
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())

#-- grid lines ---
gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,90, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

#-- set extent and coastlines----
#axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.set_extent([105,150,10,45])  # for figure
axmap.coastlines(color="k")

#-- color boundaries norm --------
#cmbnd = list(np.arange(-20,20+0.1,5))
#cmlabels = list(np.arange(-20,20+0.1,5))

cmbnd = list(np.arange(-40,-5+0.01,5)) + [-2,2] + list(np.arange(5,40+0.01,5))
cmlabels = [-40,-30,-20,-10,-5,-2,2,5,10,20,30,40]

mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

##-- contourf --------------
#X,Y = np.meshgrid(o1loncnt,o1latcnt)
##vmin, vmax = cmbnd[0], cmbnd[-1]
#im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
#print(cmbnd)

#-- pcolormesh ----
X,Y = np.meshgrid(o1lonbnd, o1latbnd)
im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm)
print(cmbnd)

#-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -1, 1) # should be adjusted considering
#plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)
plt.pcolor(X, Y, a2hatch, hatch="//", alpha=0.)
 
#-- draw colorbar ------
cax = figmap.add_axes([0.85,0.1,0.03,0.8])
cbar= plt.colorbar(im, orientation='vertical', cax=cax)
cbar.set_ticks(cmlabels)
cbar.set_ticklabels(cmlabels)

#-- vector ------------
#axmap.quiver(X, Y, a2difx, a2dify, angles="xy", units="xy", scale=1, color="gray", width=0.4)

#-- title, figure name -----
nens = a3avehis.shape[0]
figpath = figdir + '/map.difrat.speed.%04d-%04d.png'%(iY,eY)
stitle = 'TC speed diff. (%)\n' + '%04d-%04d nens:%03d'%(iY,eY,nens)

axmap.set_title(stitle)

plt.show()

plt.savefig(figpath)
print(figpath)


#****************************************************
# Figure: difference (in equivalent TC count)
#****************************************************
for scen in ["HPB","HPB_NAT"]:
    lens = list(range(1,50+1))
    a3freq = np.empty([len(lens),nyout,nxout], "float32")
    for iens,ens in enumerate(lens):
        a2count = np.zeros([nyout,nxout], "int32") 
        nstep   = 0
        for Year in lYear:
            freqbasedir = '/home/utsumi/temp/bams2020/map-freq'

            freqdir = freqbasedir + '/%s/%s-%03d'%(slabel, scen, ens)

            countpath= freqdir + '/a2count.tc.obj.%04d.npy'%(Year)
            nsteppath= freqdir + '/nstep.tc.obj.%04d.npy'%(Year)

            a2count = a2count + np.load(countpath)
            nstep   = nstep + np.load(nsteppath)

        a3freq[iens] = ma.masked_less(a2count,0) *4* 365/ float(nstep) # times / year
    if scen=="HPB":
        a2freq_his = a3freq.mean(axis=0)
    elif scen=="HPB_NAT":
        a2freq_nat = a3freq.mean(axis=0)

    if scen=="HPB":
        a2factor = 1- ma.masked_where(a2avenat==0, a2avehis)/a2avenat
        a2fig = a3freq.mean(axis=0) * a2factor  # equivalent count difference
    elif scen=="HPB_NAT":
        a2factor = ma.masked_where(a2avehis==0, a2avenat)/a2avehis - 1
        a2fig = a3freq.mean(axis=0) * a2factor  # equivalent count difference


    a2hatch = ma.masked_where(a2pv>0.05, a2fig)
    #a2hatch = ma.masked_where(a2fig<0., a2fig)
    
    #---------------------------
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
    
    #-- grid lines ---
    gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,90, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
    axmap.set_yticks(yticks, crs = ccrs.PlateCarree())
    
    #-- set extent and coastlines----
    #axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.set_extent([105,150,10,45])  # for figure
    axmap.coastlines(color="k")
    
    #-- color boundaries norm --------
    #cmbnd = list(np.arange(-1.2,1.2+0.01, 0.1))
    #cmlabels = [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2]
    
    cmbnd0 = np.array([0.1,0.2,0.3,0.5,1.0,1.5,2.0])
    cmbnd = list(-1*cmbnd0[::-1]) + [0] + list(cmbnd0)
    cmlabels = cmbnd

    mycm = 'RdBu_r'
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    
    ##-- contourf --------------
    #X,Y = np.meshgrid(o1loncnt,o1latcnt)
    ##vmin, vmax = cmbnd[0], cmbnd[-1]
    #im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
    #print(cmbnd)

    #-- pcolormesh ----
    X,Y = np.meshgrid(o1lonbnd, o1latbnd)
    im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm)
    print(cmbnd)
    
    #-- hatch --------------
    #a2hatch = ma.masked_inside(a2hatch, -1, 1) # should be adjusted considering
    #plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)
    plt.pcolor(X, Y, a2hatch, hatch="//", alpha=0.)
    
    #-- draw colorbar ------
    cax = figmap.add_axes([0.85,0.2,0.03,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)
    
    #-- vector ------------
    #axmap.quiver(X, Y, a2difx, a2dify, angles="xy", units="xy", scale=1, color="gray", width=0.4)
    
    #-- title, figure name -----
    nens = a3avehis.shape[0]
    figpath = figdir + '/map.difcount.speed.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'TC count diff. by speed change (%s-based)\n'%(scen) + '%04d-%04d nens:%03d'%(iY,eY,nens)
    
    axmap.set_title(stitle)
    
    plt.show()
    
    plt.savefig(figpath)
    print(figpath)


#****************************************************
# Figure: difference (relative to absolute TC count difference)
#****************************************************
for scen in ["HPB_NAT"]:
    a2dfreq = a2freq_his - a2freq_nat
    if scen=="HPB":
        a2factor = 1- ma.masked_where(a2avenat==0, a2avehis)/a2avenat
        a2fig = ma.masked_where(a2dfreq==0, a3freq.mean(axis=0) * a2factor) / np.abs(a2dfreq) * 100   # percent in equivalent count difference
    elif scen=="HPB_NAT":
        a2factor = ma.masked_where(a2avehis==0, a2avenat)/a2avehis - 1
        a2fig = ma.masked_where(a2dfreq==0, a3freq.mean(axis=0) * a2factor) / np.abs(a2dfreq) * 100  # percent in equivalent count difference

    a2fig = ma.masked_where(ma.masked_inside(a2dfreq, -0.5, 0.5).mask, a2fig)

    a2hatch = ma.masked_where(a2pv>0.05, a2fig)
    #a2hatch = ma.masked_where(a2fig<0., a2fig)
    
    #---------------------------
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
    axmap.set_facecolor("gray") 

    #-- grid lines ---
    gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,90, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
    axmap.set_yticks(yticks, crs = ccrs.PlateCarree())
    
    #-- set extent and coastlines----
    #axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.set_extent([105,150,10,45])  # for figure
    axmap.coastlines(color="k")
    
    #-- color boundaries norm --------
    #cmbnd = list(np.arange(-1.2,1.2+0.01, 0.1))
    #cmlabels = [-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2]

    cmbnd = np.arange(-50,50+1,10) 
    cmlabels = cmbnd

    mycm = 'RdBu_r'
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    
    #-- contourf --------------
    #X,Y = np.meshgrid(o1loncnt,o1latcnt)
    #vmin, vmax = cmbnd[0], cmbnd[-1]
    #im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")

    #-- pcolormesh ----
    X,Y = np.meshgrid(o1lonbnd, o1latbnd)
    im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm)
    print(cmbnd)
    
    #-- hatch --------------
    #plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)
    plt.pcolor(X,Y,a2hatch, hatch="//", alpha=0)
    
    #-- draw colorbar ------
    cax = figmap.add_axes([0.85,0.2,0.03,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)
    
    #-- vector ------------
    #axmap.quiver(X, Y, a2difx, a2dify, angles="xy", units="xy", scale=1, color="gray", width=0.4)
    
    #-- title, figure name -----
    nens = a3avehis.shape[0]
    figpath = figdir + '/map.difcount-frac-to-dcount.speed.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'TC %% to all dif. count by speed change (%s-based)\n'%(scen) + '%04d-%04d nens:%03d'%(iY,eY,nens)
    
    axmap.set_title(stitle)
    
    plt.show()
    
    plt.savefig(figpath)
    print(figpath)







# %%
