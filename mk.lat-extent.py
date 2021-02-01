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
#lens    = range(1,3+1)
#lens    = range(21,50+1)
lens    = range(1,50+1)

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
a1latBnd= d4PDF.LatBnd()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999
[[lllat,lllon],[urlat,urlon]] = [[0,100],[90,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

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

            a1maxlat = deque([])
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

                    amask= np.fromfile(maskpath, bool)
                    alat = np.fromfile(srcdir + "/lat.%s.bn"%(cid), "float32")
                    alon = np.fromfile(srcdir + "/lon.%s.bn"%(cid), "float32")

                    amasklat= ma.masked_outside(alat, lllat, urlat).mask
                    amasklon= ma.masked_outside(alon, lllon, urlon).mask

                    amask = amask + amasklat + amasklon

                    if (~amask).sum()==0: continue
                    if amask.sum() == 0:
                        maxlat = alat.max()
                    else:
                        maxlat = ma.masked_where(amask, alat).max()   # mask is effective for max()
    
                    a1maxlat.append(maxlat)
            a1maxlat = np.array(a1maxlat).astype("float32")

            #-- save data ---
            latbasedir = "/home/utsumi/temp/bams2020/tc-lat-extent/%s/%s.lat%03d-%03d.lon%03d-%03d"%(slabel,season,lllat,urlat,lllon,urlon)
            latdir = latbasedir + "/%s.%03d"%(scen,ens)
            util.mk_dir(latdir)
            print(latdir)
            latpath = latdir + "/lat-extent.%04d.bn"%(Year)
            a1maxlat.tofile(latpath)

#----Draw ---------------
for scen in lscen:
    adatall = deque([])
    for iens, ens in enumerate(lens):
        #-- Load --------
        latbasedir = "/home/utsumi/temp/bams2020/tc-lat-extent/%s/%s.lat%03d-%03d.lon%03d-%03d"%(slabel,season,lllat,urlat,lllon,urlon)
        latdir = latbasedir + "/%s.%03d"%(scen,ens)

        adat = np.concatenate([np.fromfile(latdir + "/lat-extent.%04d.bn"%(Year) , "float32") for Year in lYear])

        adatall.extend(adat)
        #-- histogram ---
        abnd = a1latBnd[int(nyin/2):][::4]
        acnt = (abnd[:-1] + abnd[1:])*0.5

        #----------------
        acount, _ = np.histogram(adat, abnd, density=False)
        arfreq, _ = np.histogram(adat, abnd, density=True)
        if iens==0:
            a2count = np.zeros([len(lens),len(acnt)])
            a2rfreq = np.zeros([len(lens),len(acnt)])

        a2count[iens] = acount
        a2rfreq[iens] = arfreq


    adatall = np.array(adatall)
    acountall, _ = np.histogram(adatall, abnd, density=False)
    arfreqall, _ = np.histogram(adatall, abnd, density=True)

    if scen=="HPB":
        a2count_his = a2count
        a2rfreq_his = a2rfreq
        a1countall_his = acountall
        a1rfreqall_his = arfreqall
        avehis = adatall.mean()

    elif scen=="HPB_NAT":
        a2count_nat = a2count
        a2rfreq_nat = a2rfreq
        a1countall_nat = acountall
        a1rfreqall_nat = arfreqall

        avenat = adatall.mean()

    else:
        print("check scen",scen)
        sys.exit()


#--- draw (bar histogram) -------
if figflag==True:
    fig = plt.figure(figsize=(6,6))
    for i,density in enumerate([True, False]):
        ax  = fig.add_subplot(2,1,i+1)
        if density==True:
            a1his = a1rfreqall_his
            a1nat = a1rfreqall_nat
    
        else:
            a1his = a1countall_his
            a1nat = a1countall_nat

        wbar = abnd[1] - abnd[0] 
        ax.bar(acnt, a1his, width=wbar, alpha=0.5, color="red")
        ax.bar(acnt, a1nat, width=wbar, alpha=0.5, color="blue")
    
        ax.axvline(avehis, linestyle="-", linewidth=1, color="red")
        ax.axvline(avenat, linestyle="-", linewidth=1, color="blue")
    
        ax.set_yscale("log")
        stitle = "Latitude-extent %04d-%04d ens:%03d-%03d"%(iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
    
        figpath = figdir + "/bar.lat-extent.lat%03d-%03d.lon%03d-%03d.png"%(lllat,urlat, lllon, urlon)
    
    plt.suptitle(stitle)
    plt.savefig(figpath)
    plt.show()

    print(figpath)

# %%
