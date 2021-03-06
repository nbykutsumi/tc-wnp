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
calcflag = True
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
#lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,10+1)
#lens    = range(1,50+1)
lens    = [1]

dlatsub  = 5
dlonsub  = 5
llatsub0 = np.arange(5,40+0.1,dlatsub)[::-1]
llonsub0 = np.arange(105,145+0.1,dlonsub)
lsubbox = [[[latsub0,lonsub0],[latsub0+dlatsub,lonsub0+dlonsub]] for latsub0 in llatsub0 for lonsub0 in llonsub0]

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
    a3num = np.zeros([len(lens), nyout, nxout], int32) 
    for iens, ens in enumerate(lens):
        print("load",scen,ens)
        predir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "prepos", scen, ens)
        nowdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "nowpos", scen, ens)
        nexdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "nextpos", scen, ens)

        a1pre = np.concatenate([np.load(predir + "/prepos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1now = np.concatenate([np.load(nowdir + "/nowpos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1nex = np.concatenate([np.load(nexdir + "/nextpos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])


        a1pre = ma.masked_less(a1pre, 0)
        a1nex = ma.masked_less(a1nex, 0)

        a1maskpre = a1pre.mask 
        a1masknex = a1nex.mask

        a1pre = (a1pre-1).filled(0)  # fortran index to python index
        a1nex = (a1nex-1).filled(0)  # fortran index to python index
        a1now = a1now -1

        a1ypre, a1xpre = np.unravel_index(a1pre, [nyin, nxin])
        a1ynow, a1xnow = np.unravel_index(a1now, [nyin, nxin])
        a1ynex, a1xnex = np.unravel_index(a1nex, [nyin, nxin])

        a1latpre = a1latin[a1ypre]
        a1lonpre = a1lonin[a1xpre]
        a1latnow = a1latin[a1ynow]
        a1lonnow = a1lonin[a1xnow]
        a1latnex = a1latin[a1ynex]
        a1lonnex = a1lonin[a1xnex]

        a1distpre = util.calc_dist_gc_array(a1latpre, a1latnow, a1lonpre, a1lonnow) 
        a1distnex = util.calc_dist_gc_array(a1latnex, a1latnow, a1lonnex, a1lonnow) 
        a2dist = np.array([a1distpre, a1distnex])
        a2mask = np.array([a1maskpre, a1masknex])
        a1vel = (ma.masked_where(a2mask, a2dist).mean(axis=0)) / 6.0  # km/6-hrs --> km/hr

        a1dlatpre = a1latnow - a1latpre
        a1dlonpre = a1lonnow - a1lonpre
        a1dlatnex = a1latnex - a1latnow
        a1dlonnex = a1lonnex - a1lonnow

        a1velxpre = a1vel

        a1mask = a1vel.mask
        if a1mask.sum()>0:
            a1vel   = a1vel.compressed()
            a1latnow = ma.masked_where(a1mask, a1latnow).compressed()
            a1lonnow = ma.masked_where(a1mask, a1lonnow).compressed()

        a2ave, _,_,_ = stats.binned_statistic_2d(a1latnow, a1lonnow, a1vel, statistic="mean", bins=[o1latbnd, o1lonbnd]) 
        a2sum, _,_,_ = stats.binned_statistic_2d(a1latnow, a1lonnow, a1vel, statistic="sum", bins=[o1latbnd, o1lonbnd]) 
        a2num, _,_,_ = stats.binned_statistic_2d(a1latnow, a1lonnow, a1vel, statistic="count", bins=[o1latbnd, o1lonbnd]) 
        
        a3ave[iens] = a2ave   
        a3sum[iens] = a2sum
        a3num[iens] = a2num
        plt.imshow(a2ave, origin="lower")
        plt.show()
        print(a2ave.shape)
        print(a2ave.max())
        print("lat")
        print(a1latnow)
        print(o1latbnd)
        sys.exit()

    #-- Save 3D data -----
    speeddir = speedbasedir + "/%s"%(slabel)
    util.mk_dir(speeddir)
    avepath= speeddir + "/ave.%s.%04d-%04d.npy"%(scen,iY,eY)
    sumpath= speeddir + "/sum.%s.%04d-%04d.npy"%(scen,iY,eY)
    numpath= speeddir + "/num.%s.%04d-%04d.npy"%(scen,iY,eY)

    np.save(avepath, a2ave)
    np.save(sumpath, a2sum)
    np.save(numpath, a2num)

    print(avepath)

##-- Draw ----
#a2numhis = a3numhis.sum(axis=0)
#a2numnat = a3numnam.sum(axis=0)
#a2sumhis = a3sumhis.sum(axis=0)
#a2sumnat = a3sumnat.sum(axis=0)
#a2avehis = ma.masked_where(a2numhis==0, a2sumhis) / a2numhis
#a2avenat = ma.masked_where(a2numnat==0, a2sumnat) / a2numnat



# %%
