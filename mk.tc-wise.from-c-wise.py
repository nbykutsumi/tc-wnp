# %%
#import matplotlib
#matplotlib.use('Agg')
##%matplotlib inline
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#import matplotlib.ticker as mticker

#----------------------------------
import sys, os, shutil, socket
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
from bisect import bisect_left
from collections import deque
import glob
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 1990,1990
#iY, eY = 2001,2010
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
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(21,50+1)

res     = "320x640"
noleap  = False

#lvname = ["dtlw","dtmd","dtup","lat","lon","slp","slp_mean_adj","slp_mean_box","vortlw","vortlw_max_box","wmaxlw","wmaxup","x","y"]
lvname = ["dtlw","dtmd","dtup","lat","lon","slp","slp_mean_adj","slp_mean_box","vortlw","vortlw_max_box","wmaxlw","wmaxup","x","y"]

ddtype={
"dtlw": float32,
"dtmd": float32,
"dtup": float32,
"initland": float32,
"initsst": float32,
"lat": float32,
"lon": float32,
"slp": float32,
"slp_mean_adj": float32,
"slp_mean_box": float32,
"vortlw": float32,
"vortlw_max_box": float32,
"wmaxlw": float32,
"wmaxup": float32,
"dura":int32,
"x": int32,
"y": int32,
"tcmask":bool,
}

detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))


hostname = socket.gethostname()
if hostname=="shui":
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
elif hostname=="well":
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
else:
    print("check hostname",hostname)
    sys.exit()

outbasedir = '/home/utsumi/temp/bams2020/tc-wise-WNP'
#util.mk_dir(outbaseDir)
figdir  = '/home/utsumi/temp/bams2020/fig/pdf-tc-var'
util.mk_dir(figdir)
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]

x0 = bisect_left(a1lonin, lllon)
x1 = bisect_left(a1lonin, urlon)
y0 = bisect_left(a1latin, lllat)
y1 = bisect_left(a1latin, urlat)


initx0 = bisect_left(a1lonin, 90)
initx1 = bisect_left(a1lonin, 270)
inity0 = bisect_left(a1latin, 0)
inity1 = bisect_left(a1latin, 35)


#************************************
# d4PDF (Objective detection)
#************************************
thsstdeg = 27
exrvort = 3*1.0e-5
tcrvort = 3*1.0e-5
thwcore= 0
thdura = 36
thwind = 12
thwdif = -9999



thsstk     = thsstdeg +273.15
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsstdeg*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

dhis = {}
dnat = {}
for scen in lscen:
    for ens in lens:
        for Year in lYear:
            for Mon in lMon:
                print(scen,ens,Year,Mon)
                ssearch = wsbaseDir + "/XX-%s-%03d/6hr/cwise/%04d/%02d??/initsst.*.bn"%(scen,ens,Year,Mon)
                lisstpath = sorted(glob.glob(ssearch))

                for isstpath in lisstpath:
                    fname = os.path.basename(isstpath)
                    ipos = int(fname.split(".")[-2])
                    inity,initx = np.unravel_index(ipos, [nyin, nxin])
                    if initx<initx0: continue
                    if initx>initx1: continue
                    if inity<inity0: continue
                    if inity>inity1: continue

                    idir = os.path.dirname(isstpath)
                    cid  = ".".join(fname.split(".")[-3:-1])

                    isst  = np.fromfile(isstpath, ddtype["initsst"])
                    iland = np.fromfile(idir + "/initland.%s.bn"%(cid), ddtype["initland"])

                    #-- Genesis condition ---
                    if isst  < thsstk: continue
                    if iland > 0.0: continue

                    #-- maximum wind @850hPa --
                    awmaxlw = np.fromfile(idir + "/wmaxlw.%s.bn"%(cid),ddtype["wmaxlw"])
                    nrec    = len(awmaxlw)

                    a1mask = ma.masked_less(awmaxlw,  thwind).mask

                    if a1mask.sum()==nrec: continue

                    #-- wup & wlw --
                    awmaxup = np.fromfile(idir + "/wmaxup.%s.bn"%(cid), ddtype["wmaxup"])
                    a1mask = a1mask + ma.masked_less(awmaxlw - awmaxup , thwdif).mask

                    if a1mask.sum()==nrec: continue

                    #-- thwcore  ---
                    adtlw = fromfile(idir + "/dtlw.%s.bn"%(cid), dtype=ddtype["dtlw"])
                    adtmd = fromfile(idir + "/dtmd.%s.bn"%(cid), dtype=ddtype["dtmd"])
                    adtup = fromfile(idir + "/dtup.%s.bn"%(cid), dtype=ddtype["dtup"])
                    a1mask = a1mask + ma.masked_less(adtlw+adtmd+adtup , thwcore).mask

                    if a1mask.sum()==nrec: continue
                    #-- thwovrt ----
                    avortlw = np.fromfile(idir + "/vortlw.%s.bn"%(cid), dtype=ddtype["vortlw"])
                    a1mask = a1mask + ma.masked_less(avortlw, tcrvort).mask

                    if a1mask.sum()==nrec: continue


                    #--- Save data ---
                    outdir = outbasedir + "/%s/%s.%03d/%04d%02d"%(slabel,scen,ens,Year,Mon) 
                    util.mk_dir(outdir)

                    for vname in lvname:
                        srcpath = idir   + "/%s.%s.bn"%(vname,cid)
                        outpath = outdir + "/%s.%s.bn"%(vname,cid)
                        shutil.copy(srcpath, outpath)

                    maskpath = outdir + "/tcmask.%s.bn"%(cid)
                    a1mask.astype(ddtype["tcmask"]).tofile(maskpath) 

                print(maskpath)
