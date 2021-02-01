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
from math import sin, cos, acos

#--------------------------------------
#calcflag = True
calcflag = False
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
#lens    = range(1,10+1)
lens    = range(1,50+1)
#lens    = [1]

[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]   # for loading tc-vector data

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

####****************
for scen in lscen:
    if calcflag != True: continue

    a3ave = np.zeros([len(lens), nyout, nxout], float32) 
    a3sum = np.zeros([len(lens), nyout, nxout], float32) 
    a3sumx= np.zeros([len(lens), nyout, nxout], float32) 
    a3sumy= np.zeros([len(lens), nyout, nxout], float32) 
    a3num = np.zeros([len(lens), nyout, nxout], int32) 
    for iens, ens in enumerate(lens):
        print("load",scen,ens)
        nowdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "nowpos", scen, ens)
        nexdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "nextpos", scen, ens)

        a1now = np.concatenate([np.load(nowdir + "/nowpos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1nex = np.concatenate([np.load(nexdir + "/nextpos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])

        a1nex = ma.masked_less(a1nex, 0)
        a1masknex = a1nex.mask

        a1nex = (a1nex-1).filled(0)  # fortran index to python index
        a1now = a1now -1

        a1ynow, a1xnow = np.unravel_index(a1now, [nyin, nxin])
        a1ynex, a1xnex = np.unravel_index(a1nex, [nyin, nxin])

        a1latnow = a1latin[a1ynow]
        a1lonnow = a1lonin[a1xnow]
        a1latnex = a1latin[a1ynex]
        a1lonnex = a1lonin[a1xnex]

        a1latmid = (a1latnow + a1latnex)*0.5
        a1lonmid = (a1lonnow + a1lonnex)*0.5

        a1vel  = util.calc_dist_gc_array(a1latnex, a1latnow, a1lonnex, a1lonnow)  / 6.0  # km/6hrs
        a1dlat = a1latnex - a1latnow
        a1dlon = a1lonnex - a1lonnow

        a1denom = np.sqrt( np.square(a1dlon) + np.square(a1dlat) ) 
        a1velx = a1vel * a1dlon / a1denom
        a1vely = a1vel * a1dlat / a1denom

        a1maskzero = np.logical_and((a1latnex==a1latnow),(a1lonnex==a1lonnow))

        if a1maskzero.sum()>0:
            a1vel    = ma.masked_where(a1maskzero, a1vel).filled(0)
            a1velx   = ma.masked_where(a1maskzero, a1velx).filled(0)
            a1vely   = ma.masked_where(a1maskzero, a1vely).filled(0)

        a1mask = a1masknex
        if a1mask.sum()>0:
            a1vel    = ma.masked_where(a1mask, a1vel).compressed()
            a1velx   = ma.masked_where(a1mask, a1velx).compressed()
            a1vely   = ma.masked_where(a1mask, a1vely).compressed()
            a1latmid = ma.masked_where(a1mask, a1latmid).compressed()
            a1lonmid = ma.masked_where(a1mask, a1lonmid).compressed()

        a2ave,  _,_,_ = stats.binned_statistic_2d(a1latmid, a1lonmid, a1vel, statistic="mean", bins=[o1latbnd, o1lonbnd]) 
        a2sum,  _,_,_ = stats.binned_statistic_2d(a1latmid, a1lonmid, a1vel, statistic="sum", bins=[o1latbnd, o1lonbnd]) 
        a2num,  _,_,_ = stats.binned_statistic_2d(a1latmid, a1lonmid, a1vel, statistic="count", bins=[o1latbnd, o1lonbnd]) 
        a2sumx, _,_,_ = stats.binned_statistic_2d(a1latmid, a1lonmid, a1velx, statistic="sum", bins=[o1latbnd, o1lonbnd]) 
        a2sumy, _,_,_ = stats.binned_statistic_2d(a1latmid, a1lonmid, a1vely, statistic="sum", bins=[o1latbnd, o1lonbnd]) 
        
        a3ave[iens] = a2ave   
        a3sum[iens] = a2sum
        a3sumx[iens] = a2sumx
        a3sumy[iens] = a2sumy
        a3num[iens] = a2num

    #-- Save 3D data -----
    speeddir = speedbasedir + "/%s"%(slabel)
    util.mk_dir(speeddir)
    avepath= speeddir + "/ave.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumpath= speeddir + "/sum.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumxpath= speeddir + "/sumx.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumypath= speeddir + "/sumy.%s.%04d-%04d.npy"%(scen,iY,eY)
    numpath= speeddir + "/num.%s.%04d-%04d.npy"%(scen,iY,eY)

    np.save(avepath, a3ave)
    np.save(sumpath, a3sum)
    np.save(sumxpath, a3sumx)
    np.save(sumypath, a3sumy)
    np.save(numpath, a3num)

    print(avepath)

#-- Draw ----
for scen in lscen:
    speeddir = speedbasedir + "/%s"%(slabel)
    #avepath= speeddir + "/ave.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumpath= speeddir + "/sum.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumxpath= speeddir + "/sumx.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumypath= speeddir + "/sumy.%s.%04d-%04d.npy"%(scen,iY,eY)
    numpath= speeddir + "/num.%s.%04d-%04d.npy"%(scen,iY,eY)

    #a3ave  = np.load(avepath)
    a3sum  = ma.masked_invalid(np.load(sumpath))
    a3sumx = ma.masked_invalid(np.load(sumxpath))
    a3sumy = ma.masked_invalid(np.load(sumypath))
    a3num  = ma.masked_invalid(np.load(numpath))

    a2sum = a3sum.sum(axis=0)
    a2num = a3num.sum(axis=0)
    a2ave = ma.masked_where(a2num==0, a2sum) / a2num

    a2sumx= a3sumx.sum(axis=0)
    a2sumy= a3sumy.sum(axis=0)
    a2avex= ma.masked_where(a2num==0, a2sumx) / a2num
    a2avey= ma.masked_where(a2num==0, a2sumy) / a2num

    #** Drawa ****************
    a2fig = a2ave
    a1lat = o1latcnt
    a1lon = o1loncnt
    X,Y = np.meshgrid(a1lon, a1lat)
    BBox = [[0,100],[50,150]]
    [[lllat,lllon],[urlat,urlon]] = BBox

    fig = plt.figure(figsize=(6,4))
    axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
    axmap.coastlines()

    gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks = np.arange(-180,180+1,15)
    yticks = np.arange(-90,90+1,15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
    axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

    axmap.set_extent([lllon,urlon,lllat,urlat])

    #-- Make new colormap (white at the lower end) --
    upper = matplotlib.cm.jet(np.arange(256))
    #lower = np.ones((int(256/4),4))
    lower = np.ones((int(256/8),4))
    for i in range(3):
        lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
    mycm = np.vstack(( lower, upper ))
    mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

    #-- color boundaries norm ------
    cmbnd = list(np.arange(0,40+1,2))
    #cmlabels = list(np.arange(0,50+1,10))

    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    #axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.set_extent([105,150,10,45])
    axmap.coastlines()

    #--contourf ----------
    im = axmap.contourf(X,Y, a2fig, levels=cmbnd, cmap=mycm, extend="max")
    plt.colorbar(im)

    #-- vector ------------
    axmap.quiver(X, Y, a2avex, a2avey, angles="xy", units="xy", scale=5, color="white", width=0.4)

    #-- title, figure name -----
    figpath = figdir + '/map.speed.%s.%04d-%04d.png'%(scen,iY,eY)
    stitle = 'TC speed km/h %s\n'%(scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])

    axmap.set_title(stitle)

    plt.show()

    plt.savefig(figpath)
    print(figpath)


# %%
