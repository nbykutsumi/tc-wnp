# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import scipy.stats
#----------------------------------
import sys, os, pickle
#from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import util
import calendar
#import Cyclone
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
lYear = range(iY,eY+1)
#-----------------
prj     = "d4PDF"
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(1,5+1)
#lens    = range(1,4+1)

figdir  = '/home/utsumi/temp/bams2020/fig/tc-prec'
util.mk_dir(figdir)
#----------------------------------
detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
a1latcnt  = d4PDF.Lat()
a1loncnt  = d4PDF.Lon()
a1latbnd  = d4PDF.LatBnd()
a1lonbnd  = d4PDF.LonBnd()
a1lonbnd[0] = 0
a1lonbnd[-1] = 360
miss_int= -9999
ny = len(a1latcnt)
nx = len(a1loncnt) 
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

#**************************
miss =-9999
miss_int= -9999
thsst = 27
exrvort = 3*1e-5
tcrvort = 3*1e-5
thwcore  = 0
thwind   = 12
thwdif   = -9999  # not used
thdura   = 36
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5
#************************************
# Load data
#************************************
for scen in ["HPB","HPB_NAT"]:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    a3num = np.zeros([len(lens),ny,nx], "int32")
    for iens,ens in enumerate(lens):
        print(scen,ens)
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-annmax"

        a2numTmp = np.array([np.load(precbasedir + "/%s.%03d/num.%04d.npy"%(scen,ens,Year)) for Year in lYear]).sum(axis=0)

        a3num[iens] = a2numTmp

    if scen=="HPB":
        a3dat_his = a3num
    elif scen=="HPB_NAT":
        a3dat_nat = a3num
#***********************
# Figure: ensemble mean
#***********************
a3dat_his = ma.masked_less(a3dat_his, 0)
a3dat_nat = ma.masked_less(a3dat_nat, 0)

a2dif = (a3dat_his.mean(axis=0) - a3dat_nat.mean(axis=0)) /len(lYear) *10   # times/10-year

a2tv, a2pv = scipy.stats.ttest_ind(a3dat_his.filled(np.nan), a3dat_nat.filled(np.nan), axis=0, equal_var=False, nan_policy="omit")

a2fig = a2dif
a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#a2hatch = ma.masked_where(a2fig<0., a2fig)

#-- title, figure name -----
stitle = 'diff. tc-precip-extreme counts (times/10-yr)\n' + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 
figpath = figdir + '/map.dif.tc-prec-extreme-counts.%04d-%04d.png'%(iY,eY)
#---------------------------
figmap   = plt.figure(figsize=(6,4))
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())

#-- grid lines ---
gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,90, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
axmap.set_xticks(xticks, crs=ccrs.PlateCarree())
axmap.set_yticks(yticks, crs=ccrs.PlateCarree())
#-- set extent and coastlines----
axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.coastlines(color="k")

#-- color boundaries norm --------
cmbnd = list(np.arange(-2.5,2.5+0.01,0.25))
cmlabels = list(np.arange(-2.5,2.5+0.01,1))
mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

#-- contourf --------------
X,Y = np.meshgrid(a1loncnt,a1latcnt)
#vmin, vmax = cmbnd[0], cmbnd[-1]
im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
print(cmbnd)
#-- hatch --------------
plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)

#-- draw colorbar ------
cax = figmap.add_axes([0.82,0.2,0.05,0.6])
cbar= plt.colorbar(im, orientation='vertical', cax=cax)
#cbar.set_ticks(cbar.ax.get_yticks())
cbar.set_ticks(cmlabels)
cbar.set_ticklabels(cmlabels)

#-- Tiele -----
axmap.set_title(stitle)
#-- Save ------
plt.savefig(figpath)
plt.show()
print(figpath)


###--- test map -------------
#a2his = a3rat_his.mean(axis=0).filled(0)
#a2nat = a3rat_nat.mean(axis=0).filled(0)
#a2dif = a2his - a2nat
#X,Y = np.meshgrid(a1loncnt,a1latcnt)
##clev = range(-200,200+1,20)
#clev = list(np.arange(-0.03,0.03+0.001,0.003))
#figmap = plt.figure()
#axmap  = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())
#axmap.set_extent([lllon,urlon,lllat,urlat])
#im  = plt.contourf(X,Y,a2dif, cmap="RdBu_r", levels=clev)
##im  = plt.contourf(X,Y,a2dif, cmap="RdBu_r")
#axmap.coastlines(color="k")
#plt.colorbar(im)
#plt.savefig(figdir+"/temp.dif.png")
#print("his")
#
#
#
#
#a2his = a3rat_his.mean(axis=0)*60*60*24
#a2nat = a3rat_nat.mean(axis=0)*60*60*24
#X,Y = np.meshgrid(a1loncnt,a1latcnt)
#clev = range(0,30)
#figmap = plt.figure()
#axmap  = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())
#axmap.set_extent([lllon,urlon,lllat,urlat])
#im  = plt.contourf(X,Y,a2his, clev)
##axmap.set_extent([lllon,urlon,lllat,urlat])
#axmap.coastlines(color="k")
#plt.colorbar(im)
#plt.savefig(figdir+"/temp.his.png")
#print("his")
#
#figmap = plt.figure()
#axmap  = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())
#axmap.set_extent([lllon,urlon,lllat,urlat])
#im  = plt.contourf(X,Y,a2nat, clev)
#axmap.coastlines(color="k")
#plt.colorbar(im)
#plt.savefig(figdir+"/temp.nat.png")
#print("nat")


#
## %%

# %%
