# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#----------------------------------
import sys, os, pickle
from   numpy import *
import scipy.stats
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
from bisect import bisect_left
from detect_fsub import *
import socket
#--------------------------------------
#calcflag = True
calcflag = False
figflag  = True
#iY = 1990
#eY = 2010

iY = 1990
eY = 2010

lYear = range(iY,eY+1)
lMon  = range(1,12+1)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
#lens    = list(range(16,20+1))
#lens    = list(range(22,50+1))
#lens    = list(range(36,50+1))
lens    = list(range(1,50+1))
noleap  = False
region  = "WNP"
detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'


compbasedir= '/home/utsumi/temp/bams2020/composite'
figdir  = '/home/utsumi/temp/bams2020/fig/map-prec'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640
miss =-9999
miss_int= -9999

#************************************
for scen in lscen:
    precbasedir  = "/home/utsumi/temp/bams2020/tc-prec-all"
    bboxpath = precbasedir + "/HPB.001/bbox.%s.npy"%(region)
    yxbboxpath = precbasedir + "/HPB.001/yxbbox.%s.npy"%(region)
    bbox   = np.load(bboxpath)
    yxbbox = np.load(yxbboxpath)

    [[lllat,lllon],[urlat,urlon]] = bbox
    [[y0,x0],[y1,x1]] = yxbbox
    nyreg = y1-y0+1
    nxreg = x1-x0+1

    a3max = np.empty([len(lens), nyreg,nxreg], "float64")
    for iens, ens in enumerate(lens):
        print(scen,ens)
        a3tmp = np.empty([len(lYear),nyreg,nxreg], "float64")
        for itmp,Year in enumerate(lYear):
            precpath = precbasedir + "/%s.%03d/prec-tc.%s.%03d.npy"%(scen,ens,region,Year)
        
            a3prec = np.load(precpath)

            a3tmp[itmp]= a3prec.max(axis=0)

        a3max[iens] = a3tmp.max(axis=0)

    if scen == "HPB":
        a3dat_his = a3max

    elif scen=="HPB_NAT":
        a3dat_nat = a3max

#***********************
# Figure: ensemble mean
#***********************
a3dat_his = ma.masked_less(a3dat_his, 0)
a3dat_nat = ma.masked_less(a3dat_nat, 0)

print(a3dat_his.max())

a2fig = (np.median(a3dat_his, axis=0) - np.median(a3dat_nat, axis=0))*60*60*24   # mm/day

a2tv, a2pv = scipy.stats.ttest_ind(a3dat_his.filled(np.nan), a3dat_nat.filled(np.nan), axis=0, equal_var=False, nan_policy="omit")

a2hatch = ma.masked_where(a2pv>0.05, a2fig)
#a2hatch = ma.masked_where(a2fig<0., a2fig)

#-- title, figure name -----
figpath = figdir + '/map.dif.prec-tc-max.%s.%04d-%04d.png'%(scen,iY,eY)
stitle = 'dif. median tc max prec. mm/day\n' + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1]) 

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
dbboxfig = {"WNP":[[0,100],[50,150]]}
[[lllatfig,lllonfig],[urlatfig,urlonfig]] = dbboxfig[region]

axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
axmap.coastlines(color="k")

#-- color boundaries norm --------
cmbnd = list(range(-100,100+1,10))
cmlabels = list(range(-100,100+1,20))
mycm = 'RdBu_r'

cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = matplotlib.colors.ListedColormap(cmaplist)
norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

#-- contourf --------------
a1latcnt = a1lat[y0:y1+1]
a1loncnt = a1lon[x0:x1+1]
print(y0,y1)
print(lllat,urlat)
print(a1latcnt)
X,Y = np.meshgrid(a1loncnt,a1latcnt)
print(Y)
#vmin, vmax = cmbnd[0], cmbnd[-1]
im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")
print(cmbnd)
#-- hatch --------------
#a2hatch = ma.masked_inside(a2hatch, -1, 1) # should be adjusted considering 
plt.contourf(X, Y, a2hatch, hatches=["//"], alpha=0.)

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




# %%
