## %%
#import matplotlib
#matplotlib.use('Agg')
#%matplotlib inline
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
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])
season ="ALL"
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,7+1)
lens    = [1]
#lens    = range(21,50+1)
#lens    = range(1,50+1)
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
    for ens in lens:
        for Year in lYear:
            print(scen,ens,Year)

            ldat = deque([])
            llmlat = deque([])
            llmlon = deque([])
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
                    else:
                        avartmp = ma.masked_where(amask, avar).compressed()

                    imax = np.argmax(avartmp)
                    ldat.append(avartmp[imax])
                    llmlat.append(alat[imax])
                    llmlon.append(alon[imax])

            ldat = np.array(ldat)
            llmlat=np.array(llmlat)
            llmlon=np.array(llmlon)
            #-- save --
            lmbasedir = '/home/utsumi/temp/bams2020/lm-WNP'
            lmdir = lmbasedir + "/%s/%s.lat%03d-%03d.lon%03d-%03d/%s.%03d/"%(slabel,season,lllat,urlat,lllon,urlon,scen,ens) 
            lmpath    = lmdir + "/lm.%s.%04d.bn"%(vname,Year)
            latpath = lmdir + "/lat.%s.%04d.bn"%(vname,Year)
            lonpath = lmdir + "/lon.%s.%04d.bn"%(vname,Year)
            util.mk_dir(lmdir)
            ldat.astype(ddtype[vname]).tofile(lmpath)  
            llmlat.astype(ddtype["lat"]).tofile(latpath) 
            llmlon.astype(ddtype["lon"]).tofile(lonpath)
            print(lonpath)
#        #-- histogram ---
#        if vname=="dslp":
#            abnd = np.arange(0,25,0.5)
#            acnt = (abnd[:-1] + abnd[1:])*0.5
#
#        elif vname=="wmaxlw":
#            abnd = np.arange(12,100,2)
#            acnt= (abnd[:-1] + abnd[1:])*0.5
#        else:
#            print("check vname",vname)
#            sys.exit()
#        #----------------
#        acount, _ = np.histogram(ldat, abnd, density=False)
#        arfreq, _ = np.histogram(ldat, abnd, density=True)
#
#        dcount[ens] = acount
#        drfreq[ens] = arfreq
#
#    ldatall = np.array(ldatall)
#    if scen=="HPB":
#        acount_his = dcount
#        arfreq_his = drfreq
#        avehis = ldatall.mean()
#
#    elif scen=="HPB_NAT":
#        acount_nat = dcount
#        arfreq_nat = drfreq
#        avenat = ldatall.mean()
#
#    else:
#        print("check scen",scen)
#        sys.exit()
#
#
##--- draw -------
#fig = plt.figure(figsize=(6,6))
#for i,density in enumerate([True, False]):
#    ax  = fig.add_subplot(2,1,i+1)
#    if density==True:
#        a2his = a2rfreq_his
#        a2nat = a2rfreq_nat
#
#    else:
#        a2his = a2count_his
#        a2nat = a2count_nat
#
#    a1his = a2his.mean(axis=0)
#    a1nat = a2nat.mean(axis=0)
#
#    a1hisup = np.percentile(a2his, 95, axis=0)
#    a1hislw = np.percentile(a2his, 5, axis=0)
#    a1natup = np.percentile(a2nat, 95, axis=0)
#    a1natlw = np.percentile(a2nat, 5, axis=0)
#
#    ax.plot(cmcnt, a1his, "-",linewidth=2,color="red")
#    ax.plot(cmcnt, a1nat, "-",linewidth=2,color="blue")
#
#    ax.plot(cmcnt, a1hisup, "-", linewidth=0.5, color="red")
#    ax.plot(cmcnt, a1hislw, "-", linewidth=0.5, color="red")
#    ax.plot(cmcnt, a1natup, "-", linewidth=0.5, color="blue")
#    ax.plot(cmcnt, a1natlw, "-", linewidth=0.5, color="blue")
#
#    ax.axvline(avehis, linestyle="-", linewidth=1, color="red")
#    ax.axvline(avenat, linestyle="-", linewidth=1, color="blue")
#
#    ax.set_yscale("log")
#    stitle = "%s %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
#
#    figpath = figdir + "/pdf.lm.%s.lat%03d-%03d.lon%03d-%03d.png"%(vname,lllat,urlat, lllon, urlon)
#
#plt.suptitle(stitle)
#plt.savefig(figpath)
#plt.show()



# %%
