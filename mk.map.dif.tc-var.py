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
model   = "__"
expr    = 'XX'
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(1,5+1)

#tcvar = "dslp"
tcvar = "wmaxlw"

wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outbaseDir = '/home/utsumi/temp/bams2020/map-tc-%s'%(tcvar)
util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/map-tc-var'
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
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

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

    a3sum = np.empty([len(lens), ny,nx], "float64")
    a3num = np.empty([len(lens), ny,nx], "int32")
    for iens, ens in enumerate(lens):
        print(scen,ens)
        outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)

        a3sum[iens] = np.array([np.load(outDir + "/a2sum.tc.obj.%04d.npy"%(Year)) for Year in lYear]).sum(axis=0)
        a3num[iens] = np.array([np.load(outDir + "/a2num.tc.obj.%04d.npy"%(Year)) for Year in lYear]).sum(axis=0)

    a3ave = (ma.masked_where(a3num==0, a3sum) / a3num).filled(-9999.)

    if tcvar=="dslp":
        a3ave = a3ave * 0.01  # hPa

    if scen=="HPB":
        a3dat_his = a3ave
    elif scen=="HPB_NAT":
        a3dat_nat = a3ave

#***********************
# Figure: ensemble mean
#***********************
a2dif = ma.masked_less(a3dat_his,0).filled(0.0).mean(axis=0) - ma.masked_less(a3dat_nat,0).filled(0.0).mean(axis=0)

a2tv, a2pv = scipy.stats.ttest_ind(a3dat_his, a3dat_nat, axis=0, equal_var=False, nan_policy="omit")

a2fig = a2dif
a2hatch = ma.masked_where(a2pv>0.05, a2fig)

#-- title, figure name -----
if tcvar=="dslp":
    svar="(hPa)"
elif tcvar=="wmaxlw":
    svar="(m/s)"
else:
    print("check tcvar",tcvar)
    sys.exit()

stitle = 'diff. (HIST - NAT) %s\n%04d-%04d'%(svar,iY,eY) 
figpath = figdir + '/map.dif.tc-%s.%04d-%04d.png'%(tcvar,iY,eY)
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
if tcvar=="dslp":
    cmbnd = list(np.arange(-1,1+0.01, 0.1))
    cmlabels = list(np.arange(-1,1+0.01, 0.2))
elif tcvar=="wmaxlw":
    cmbnd = list(np.arange(-6,6+0.01, 1))
    cmlabels = list(np.arange(-6,6+0.01, 1))


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
plt.pcolor(X, Y, a2hatch, hatch="/", alpha=0.)

#-- draw colorbar ------
cax = figmap.add_axes([0.82,0.2,0.05,0.6])
cbar= plt.colorbar(im, orientation='vertical', cax=cax)
cbar.set_ticks(cmlabels)
cbar.set_ticklabels(cmlabels)

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


# %%
