# %%
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
%matplotlib inline
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

import numpy as np
import scipy
import d4PDF
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from importlib import import_module
from numpy import ma
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen   = ['HPB','HPB_NAT'] # run={expr}-{scen}-{ens}

#lens    = list(range(1,50+1))
lens    = list(range(1,8+1))
#lens    = [1]
miss_out = -9999.

iY,eY = 1990,2010
lYear  = range(iY,eY+1)
ny,nx = 320, 640

BBox = [[0,100],[50,150]]
[[lllat,lllon],[urlat,urlon]] = BBox


hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'

detectName = 'wsd_d4pdf_20201209-py38'

a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)

#--- load -----
maxbasedir = "/home/utsumi/temp/bams2020/max-prec-d"

scen = "HPB"
a3max = ma.asarray([np.load(maxbasedir + "/%s.%03d/ann-max-prec.%04d.npy"%(scen, ens, Year)) for Year in lYear for ens in lens])
a2his = ma.masked_less(a3max,0).mean(axis=0)

scen = "HPB_NAT"
a3max = ma.asarray([np.load(maxbasedir + "/%s.%03d/ann-max-prec.%04d.npy"%(scen, ens, Year)) for Year in lYear for ens in lens])
a2nat = ma.masked_less(a3max,0).mean(axis=0)

a2dif = a2his - a2nat
a2dif = a2dif * 60*60*24.

print(a2his.shape)
#*** Draw map *************
a1lat = d4PDF.Lat()
a1lon = d4PDF.Lon()
X,Y = np.meshgrid(a1lon, a1lat)

fig = plt.figure(figsize=(6,4))
ax  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks = np.arange(-180,180+1,15)
yticks = np.arange(-90,90+1,15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
ax.set_xticks(xticks, crs = ccrs.PlateCarree())
ax.set_yticks(yticks, crs = ccrs.PlateCarree())

ax.set_extent([lllon,urlon,lllat,urlat])
#clevs = range(-20,20+1,4)
clevs = range(-100,100+1,10)
im = ax.contourf(X,Y, a2dif, clevs, cmap="bwr")
plt.colorbar(im)
plt.show()


# %%
