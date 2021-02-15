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
import myfunc.util as util
import calendar
#--------------------------------------
lYear0 = range(1961,1987+1)
#lYear1 = range(1988,2010+1)

#lYear0 = range(1971,1990+1)
#lYear1 = range(1991,2010+1)
#lYear0 = range(1980,1990+1)
#lYear1 = range(1991,2010+1)
#lYear0 = range(1961,1985+1)
#lYear1 = range(1991,2015+1)
#lYear0 = range(1980,1997+1)
#lYear1 = range(1998,2015+1)
#lYear0 = range(1961,1985+1)
#lYear1 = range(1986,2010+1)
#lYear0 = range(1951,1970+1)
#lYear1 = range(1996,2015+1)
#lYear0 = range(1951,1970+1)
#lYear1 = range(1986,2015+1)
#lYear0 = range(1951,1970+1)
#lYear1 = range(1971,1990+1)




iY0,eY0 = lYear0[0],lYear0[-1]
iY1,eY1 = lYear1[0],lYear1[-1]

#-----------------
outbaseDir = '/home/utsumi/temp/bams2020/map-freq'
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-freq'
util.mk_dir(figdir)
#----------------------------------
miss_int= -9999
dgrid = 5
a1latbnd  = np.arange(-90,90+0.01, dgrid)
a1lonbnd  = np.arange(0,360+0.01, dgrid)
ny = len(a1latbnd) - 1
nx = len(a1lonbnd) - 1
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]
[[lllat,lllon],[urlat,urlon]] = [[10,105],[45,150]]

#**************************

#************************************
# d4PDF (Objective detection)
#************************************
freqbasedir = '/home/utsumi/temp/bams2020/map-freq'
freqdir = freqbasedir + "/bst"

a3count0 = np.array([ np.load(freqdir + '/a2count.tc.bst.%04d.%02d.npy'%(Year,Mon)) for Mon in range(1,12+1) for Year in lYear0])

nstep0 = np.array([np.load(freqdir + '/nstep.tc.bst.%04d.%02d.npy'%(Year,Mon)) for Mon in range(1,12+1) for Year in lYear0]).sum()

a3count1 = np.array([ np.load(freqdir + '/a2count.tc.bst.%04d.%02d.npy'%(Year,Mon)) for Mon in range(1,12+1) for Year in lYear1])

nstep1 = np.array([np.load(freqdir + '/nstep.tc.bst.%04d.%02d.npy'%(Year,Mon)) for Mon in range(1,12+1) for Year in lYear1]).sum()

a2freq0 = a3count0.sum(axis=0) / nstep0*365*4
a2freq1 = a3count1.sum(axis=0) / nstep1*365*4

#***********************
# Figure: ensemble mean (count)
#***********************
a2dif = a2freq1 - a2freq0

a2tv, a2pv = scipy.stats.ttest_ind(a3count1/nstep1, a3count0/nstep0, axis=0, equal_var=False, nan_policy="omit")

a2fig = a2dif
#a2hatch = ma.masked_where(a2pv>0.05, a2fig)
a2hatch = ma.masked_where(a2pv>0.1, a2fig)

#-- title, figure name -----
stitle = 'diff. count/year %04d-%04d to %04d-%04d'%(iY0,eY0,iY1,eY1) 
figpath = figdir + '/map.dif.freq.tc.bst.%04d-%04d.to.%04d-%04d.png'%(iY0,eY0,iY1,eY1)
#---------------------------
figmap   = plt.figure(figsize=(6,4))
axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())

#-- grid lines ---
gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,90, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
#-- set extent and coastlines----
axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.coastlines(color="k")

#-- color boundaries norm --------
#cmbnd = list(range(-45,45+1,10))
cmbnd = list(np.arange(-4.5,4.5+0.1,1))
mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

#-- pcolormesh --------------
X,Y = np.meshgrid(a1lonbnd,a1latbnd)
vmin, vmax = cmbnd[0], cmbnd[-1]
im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, vmin=vmin,vmax=vmax, norm=norm)

#-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -5, 5) # should be adjusted considering 
a2hatch = ma.masked_inside(a2hatch, -0.5, 0.5) # should be adjusted considering 
plt.pcolor(X, Y, a2hatch, hatch="/", alpha=0.)

#-- draw colorbar ------
cax = figmap.add_axes([0.79,0.2,0.05,0.6])
cbar= plt.colorbar(im, orientation='vertical', cax=cax)
extend = "both"
if cmbnd is not None:
    cbar.set_ticks(cbar.ax.get_yticks())
    if extend=='both':
        cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:-1] + [""])

    if extend=='min':
        cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:])

    if extend=='max':
        cbar.set_ticklabels(list(cbar.ax.get_yticklabels())[:-1]+ [""])

    else:
        pass

#-- Tiele -----
axmap.set_title(stitle)
#-- Save ------
plt.savefig(figpath)
plt.show()
print(figpath)


##***********************
## Figure: ensemble mean (fractional change (%))
##***********************
#a2difrat = ma.masked_where(a3freq_nat.mean(axis=0)==0, a3freq_his.mean(axis=0) - a3freq_nat.mean(axis=0)) / a3freq_nat.mean(axis=0) * 100 # (%)
#
#a2difrat = ma.masked_where(ma.masked_inside(a2dif,-0.5,0.5).mask, a2difrat)
#
#a2tv, a2pv = scipy.stats.ttest_ind(a3freq_his, a3freq_nat, axis=0, equal_var=False, nan_policy="omit")
#
#a2fig = a2difrat
#a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#
##-- title, figure name -----
#stitle = 'diff. (HIST - NAT) (%) ' + '%04d-%04d'%(iY,eY) 
#figpath = figdir + '/map.difrat.freq.tc.obj.%04d-%04d.png'%(iY,eY)
##---------------------------
#figmap   = plt.figure(figsize=(6,4))
#axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())
#axmap.set_facecolor("gray")  # set background color
##-- grid lines ---
#gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
#xticks   = np.arange(-180, 180+1, 15)
#yticks   = np.arange(-90,90, 15)
#gl.xlocator = mticker.FixedLocator(xticks)
#gl.ylocator = mticker.FixedLocator(yticks)
##-- set extent and coastlines----
#axmap.set_extent([lllon,urlon,lllat,urlat])
#axmap.coastlines(color="k")
#
##-- color boundaries norm --------
#cmbnd = list(np.arange(-40,40+1,5))
#mycm = 'RdBu_r'
#
#cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
#cmaplist = [cmap(i) for i in range(cmap.N)]
#cmap = matplotlib.colors.ListedColormap(cmaplist)
#norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
#
##-- pcolormesh --------------
#X,Y = np.meshgrid(a1lonbnd,a1latbnd)
#vmin, vmax = cmbnd[0], cmbnd[-1]
#im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, vmin=vmin,vmax=vmax, norm=norm)
#
##-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -5, 5) # should be adjusted considering 
#plt.pcolor(X, Y, a2hatch, hatch="/", alpha=0.)
#
##-- draw colorbar ------
#cax = figmap.add_axes([0.79,0.2,0.05,0.6])
#cbar= plt.colorbar(im, orientation='vertical', cax=cax)
#extend = "both"
#if cmbnd is not None:
#    cbar.set_ticks(cbar.ax.get_yticks())
#    if extend=='both':
#        cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:-1] + [""])
#
#    if extend=='min':
#        cbar.set_ticklabels([""] + list(cbar.ax.get_yticklabels())[1:])
#
#    if extend=='max':
#        cbar.set_ticklabels(list(cbar.ax.get_yticklabels())[:-1]+ [""])
#
#    else:
#        pass
#
##-- Tiele -----
#axmap.set_title(stitle)
##-- Save ------
#plt.savefig(figpath)
#plt.show()
#print(figpath)
#
#
#
#
#
#
##
### %%
#
## %%

# %%
