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
figmean = True
#figmean = False

#iY, eY = 1980,2010
iY, eY = 1990,2010
lYear = range(iY,eY+1)
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(1,5+1)
#lens    = range(3,9+1)
outbaseDir = '/home/utsumi/temp/bams2020/map-genesis'
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-genesis'

#----------------------------------
miss_int= -9999
dgrid = 5
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]  # data
a1latbnd  = np.arange(lllat,urlat+0.01, dgrid)
a1lonbnd  = np.arange(lllon,urlon+0.01, dgrid)
nysub = len(a1latbnd) - 1
nxsub = len(a1lonbnd) - 1
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]
[[lllatfig,lllonfig],[urlatfig,urlonfig]] = [[10,105],[45,150]]

#**************************

#************************************
# d4PDF (Objective detection)
#************************************
thsst  = 27
exrvort = 3*1.0e-5
tcrvort = 3*1.0e-5
thwcore= 0
thdura = 36
thwind = 12
thwdif = -9999

exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

for scen in ["HPB","HPB_NAT"]:
    a3freq = np.zeros([len(lens),nysub,nxsub],float32)

    for iens,ens in enumerate(lens):
        outdir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
        a2freqtmp = np.array([np.load(outdir + "/num.%04d.npy"%(Year)).astype("int32") for Year in lYear]).sum(axis=0)
        a2freqtmp = a2freqtmp.astype("float32") / len(lYear) # count per year

        a3freq[iens] = a2freqtmp


    if scen=="HPB":
        a3freq_his = a3freq
    elif scen=="HPB_NAT":
        a3freq_nat = a3freq

#***********************
# Figure: ensemble mean (count)
#***********************
a2dif = a3freq_his.mean(axis=0) - a3freq_nat.mean(axis=0)

a2tv, a2pv = scipy.stats.ttest_ind(a3freq_his, a3freq_nat, axis=0, equal_var=False, nan_policy="omit")

a2fig = a2dif
a2hatch = ma.masked_where(a2pv>0.05, a2fig)

#-- title, figure name -----
stitle = 'diff. (HIST - NAT) count/year ' + '%04d-%04d'%(iY,eY) 
figpath = figdir + '/map.dif.genesis.obj.png'
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
axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
axmap.coastlines(color="k")

#-- color boundaries norm --------
#cmbnd = list(range(-45,45+1,10))
cmbnd = list(np.arange(-1.1,1.1+0.01,0.2))
cmlabels = [-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9]

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

cbar.set_ticks(cmlabels)
cbar.set_ticklabels(cmlabels)

#-- Tiele -----
axmap.set_title(stitle)
#-- Save ------
plt.savefig(figpath)
plt.show()
print(figpath)



# %%
