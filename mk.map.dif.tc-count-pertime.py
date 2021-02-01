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
from bisect import bisect_left
#import Cyclone
#--------------------------------------
#calcflag = True
calcflag = False
figflag  = True
radkm = 500 # km
wday  = 1  # days
#wday  = 3  # days
region= "WNP"
#iY, eY = 1980,2010
iY, eY = 1990,2010
lYear = range(iY,eY+1)
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(1,3+1)
#lens    = range(3,9+1)

figdir  = '/home/utsumi/temp/bams2020/fig/map-tc-count-perday'
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
#print(a1latbnd)
#print(len(a1latbnd))
#print("")
#print("")
#print(a1lonbnd)
#print(len(a1lonbnd))
#print("")
#sys.exit()
miss_int= -9999
ny = len(a1latcnt)
nx = len(a1loncnt) 
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]
dbbox = {"WNP":[[0,100],[50,150]]}

[[lllat,lllon],[urlat,urlon]] = dbbox[region]

y0 = bisect_left(a1latcnt, lllat)
y1 = bisect_left(a1latcnt, urlat)
x0 = bisect_left(a1loncnt, lllon)
x1 = bisect_left(a1loncnt, urlon)
nyreg = y1-y0+1
nxreg = x1-x0+1
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

#***********************
# Figure: ensemble mean
#***********************
if figflag != True: sys.exit()

for scen in ["HPB","HPB_NAT"]:
    a3ave = np.full([len(lens),nyreg,nxreg], -9999).astype("float64")

    for iens,ens in enumerate(lens):
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-%04dkm"%(radkm)
        precdir = precbasedir + "/%s.%03d"%(scen,ens)

        a3pass = np.full([len(lYear),nyreg,nxreg], -9999)
        a3even = np.full([len(lYear),nyreg,nxreg], -9999)

        for i,Year in enumerate(lYear):
            print(ens,Year)
            a3count = ma.masked_equal(np.load(precdir + "/timesper-%dday.%04d.npy"%(wday,Year)), 0)

            a2passTmp = a3count.sum(axis=0)
            a2evenTmp = a3count.count(axis=0)

            a3pass[i] = a2passTmp
            a3even[i] = a2evenTmp

        a2pass   = a3pass.sum(axis=0)
        a2even   = a3even.sum(axis=0)
        a2aveTmp = (ma.masked_where(a2even==0, a2pass.astype("float64")) / a2even).filled(np.nan)


        a3ave[iens] = a2aveTmp

    if scen =="HPB":
        a3ave_his = ma.masked_invalid(a3ave, 0)
    elif scen=="HPB_NAT":
        a3ave_nat = ma.masked_less(a3ave, 0)

a2dif = a3ave_his.mean(axis=0) - a3ave_nat.mean(axis=0)  # times per time

a2tv, a2pv = scipy.stats.ttest_ind(a3ave_his.filled(np.nan), a3ave_nat.filled(np.nan), axis=0, equal_var=False, nan_policy="omit")

a2fig = a2dif
a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#a2hatch = ma.masked_where(a2fig<0., a2fig)

#-- title, figure name -----
stitle = 'dif count per-%dday (HIST - NAT)\n'%(wday) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 
figpath = figdir + '/map.dif.tc-count-per%dday.obj.%04d-%04d.png'%(wday,iY,eY)
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
axmap.set_extent([lllon,urlon,lllat,urlat])
axmap.coastlines(color="k")

#-- color boundaries norm --------
cmbnd = list(np.arange(-0.8,0.8+0.01,0.1))
cmlabels = [-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]
mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

#-- contourf --------------
X,Y = np.meshgrid(a1loncnt[x0:x1+1],a1latcnt[y0:y1+1])
#vmin, vmax = cmbnd[0], cmbnd[-1]
im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
print(cmbnd)
#-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -1, 1) # should be adjusted considering 
plt.contourf(X, Y, a2hatch, hatches=["///"], alpha=0.)

#-- draw colorbar ------
cax = figmap.add_axes([0.82,0.2,0.05,0.6])
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
