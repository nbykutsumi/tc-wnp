# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

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
from scipy import stats
#--------------------------------------
#calcflag = True
calcflag = False
figflag = True
#figflag = False

#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])
season ="ALL"
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = [1,2,3]
#lens    = range(21,50+1)
lens    = range(1,50+1)
vname = "wmaxlw"
#vname = "dslp"

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
"dslp":float32,
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
[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

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

for scen in lscen:
    if calcflag !=True: continue

    for iens,ens in enumerate(lens):
        for Year in lYear:
            print(scen,ens,Year)
            a1dat_all = deque([])
            a1lat_all = deque([])

            for Mon in lMon:
                ibasedir = '/home/utsumi/temp/bams2020/tc-wise-WNP'
                idir = ibasedir + "/%s/%s.%03d/%04d%02d"%(slabel,scen,ens,Year,Mon) 
                ssearch = idir + "/tcmask.*.bn"
                lmaskpath = sorted(glob.glob(ssearch))
                for maskpath in lmaskpath:
                    #print(maskpath)
                    srcdir= os.path.dirname(maskpath)
                    fname = os.path.basename(maskpath)
                    ipos = int(fname.split(".")[-2])

                    cid  = ".".join(fname.split(".")[-3:-1])

                    amask= np.fromfile(maskpath, ddtype["tcmask"])
                    alat = np.fromfile(srcdir + "/lat.%s.bn"%(cid), ddtype["lat"])
                    alon = np.fromfile(srcdir + "/lon.%s.bn"%(cid), ddtype["lon"])

                    if vname=="dslp":
                        aslp     = np.fromfile(srcdir + "/slp.%s.bn"%(cid), ddtype["slp"])
                        aslp_adj = np.fromfile(srcdir + "/slp_mean_adj.%s.bn"%(cid), ddtype["slp"])
                        avar = aslp_adj - aslp

                    else:
                        avar = np.fromfile(srcdir + "/%s.%s.bn"%(vname,cid), ddtype[vname])

                    amasklat= ma.masked_outside(alat, lllat, urlat).mask
                    amasklon= ma.masked_outside(alon, lllon, urlon).mask

                    amask = amask + amasklat + amasklon

                    if (~amask).sum()==0: continue
                    if amask.sum() == 0:
                        avartmp = avar
                        alattmp = alat

                    else:
                        avartmp = ma.masked_where(amask, avar).compressed()
                        alattmp = ma.masked_where(amask, alat).compressed()

                    a1dat_all.extend(avartmp)
                    a1lat_all.extend(alattmp)

            a1dat_all = np.array(a1dat_all)
            a1lat_all = np.array(a1lat_all)

            a1latbnd = d4PDF.LatBnd()
            a1latcnt = d4PDF.Lat()
            a1sum, _, _ = stats.binned_statistic(x=a1lat_all, values=a1dat_all, statistic="sum", bins=a1latbnd)
            a1num, _, _ = stats.binned_statistic(x=a1lat_all, values=a1dat_all, statistic="count", bins=a1latbnd)

            ##-- save data ---
            datbasedir = "/home/utsumi/temp/bams2020/tc-lat-var/%s/%s.lat%03d-%03d.lon%03d-%03d"%(slabel,season,lllat,urlat,lllon,urlon)
            datdir = datbasedir + "/%s.%03d"%(scen,ens)
            util.mk_dir(datdir)
            sumpath = datdir + "/sum.%s.%04d.bn"%(vname,Year)
            numpath = datdir + "/num.%s.%04d.bn"%(vname,Year)
            a1sum.astype("float32").tofile(sumpath)
            a1num.astype("int32").tofile(numpath) 
            print(datdir)

##---- plot ---------------
for scen in lscen:
    if figflag !=True: continue

    a2ave = np.zeros([len(lens), nyin],float64)
    a2num = np.zeros([len(lens), nyin],int32)

    for iens, ens in enumerate(lens):
        #-- Load --------
        datbasedir = "/home/utsumi/temp/bams2020/tc-lat-var/%s/%s.lat%03d-%03d.lon%03d-%03d"%(slabel,season,lllat,urlat,lllon,urlon)
        datdir = datbasedir + "/%s.%03d"%(scen,ens)

        a2sumTmp = np.array([np.fromfile(datdir + "/sum.%s.%04d.bn"%(vname,Year) , "float32").astype("float64") for Year in lYear])

        a2numTmp = np.array([np.fromfile(datdir + "/num.%s.%04d.bn"%(vname,Year) , "int32") for Year in lYear])

        print(a2numTmp.shape)
        a1sum = a2sumTmp.sum(axis=0)
        a1num = a2numTmp.sum(axis=0)

        a1ave = (ma.masked_where(a1num==0, a1sum) / a1num).filled(-9999.)
        a2ave[iens] = a1ave
        a2num[iens] = a1num

    if scen=="HPB":
        a2ave_his = a2ave
        a2num_his = a2num

    elif scen=="HPB_NAT":
        a2ave_nat = a2ave
        a2num_nat = a2num

    else:
        print("check scen",scen)
        sys.exit()

#--- draw -------
if figflag==True:
    fig = plt.figure(figsize=(6,6))
    for i,vartype in enumerate(["ave","num"]):
        ax  = fig.add_subplot(2,1,i+1)
        if vartype=="ave":
            a2his = a2ave_his
            a2nat = a2ave_nat
    
        else:
            a2his = a2num_his
            a2nat = a2num_nat

        a2his = ma.masked_equal(a2his, -9999.)
        a2nat = ma.masked_equal(a2nat, -9999.)

        a1his = a2his.mean(axis=0)
        a1nat = a2nat.mean(axis=0)

        acnt = d4PDF.Lat() 
    
        a1hisup = np.nanpercentile(a2his.filled(np.nan), 95, axis=0)
        a1hislw = np.nanpercentile(a2his.filled(np.nan), 5, axis=0)
        a1natup = np.nanpercentile(a2nat.filled(np.nan), 95, axis=0)
        a1natlw = np.nanpercentile(a2nat.filled(np.nan), 5, axis=0)
    
        ax.plot(acnt, a1his, "-",linewidth=2,color="red")
        ax.plot(acnt, a1nat, "-",linewidth=2,color="blue")
    
        ax.plot(acnt, a1hisup, "-", linewidth=0.5, color="red")
        ax.plot(acnt, a1hislw, "-", linewidth=0.5, color="red")
        ax.plot(acnt, a1natup, "-", linewidth=0.5, color="blue")
        ax.plot(acnt, a1natlw, "-", linewidth=0.5, color="blue")

        ax.set_xlim([0,50])    
        stitle = "%s %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
    
        figpath = figdir + "/lat-%s.lat%03d-%03d.lon%03d-%03d.png"%(vname,lllat,urlat, lllon, urlon)
    
    plt.suptitle(stitle)
    plt.savefig(figpath)
    plt.show()



# %%
