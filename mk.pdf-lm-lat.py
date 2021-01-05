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
#---- PDF ---------------
for scen in lscen:
    adatall = deque([])
    for iens, ens in enumerate(lens):
        #-- Load --------
        lmbasedir = "/home/utsumi/temp/bams2020/tc-lm/%s/%s.lat%03d-%03d.lon%03d-%03d"%(slabel,season,lllat,urlat,lllon,urlon)
        lmdir = lmbasedir + "/%s.%03d"%(scen,ens)

        adat = np.concatenate([np.fromfile(lmdir + "/lat.%s.%04d.bn"%(vname,Year) , ddtype[vname]) for Year in lYear])

        adatall.extend(adat)
        #-- histogram ---
        abnd = d4PDF.LatBnd()[160:280]
        acnt = (abnd[:-1] + abnd[1:])*0.5

        #----------------
        acount, _ = np.histogram(adat, abnd, density=False)
        arfreq, _ = np.histogram(adat, abnd, density=True)

        if iens==0:
            a2count = np.zeros([len(lens),len(acnt)])
            a2rfreq = np.zeros([len(lens),len(acnt)])

        a2count[iens] = acount
        a2rfreq[iens] = arfreq

    print(scen,a2count.mean())

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

print(a2count_his.mean(), a2count_nat.mean())


#--- draw (bar) -------
fig = plt.figure(figsize=(6,6))
for i,density in enumerate([True, False]):
    ax  = fig.add_subplot(2,1,i+1)
    if density==True:
        a1his = a1rfreqall_his
        a1nat = a1rfreqall_nat

    else:
        a1his = a1countall_his
        a1nat = a1countall_nat

    ax.bar(acnt, a1his, width=0.6, alpha=0.3, color="red")
    ax.bar(acnt, a1nat, width=0.6, alpha=0.3, color="blue")

    ax.axvline(avehis, linestyle="-", linewidth=1, color="red")
    ax.axvline(avenat, linestyle="-", linewidth=1, color="blue")

    ax.set_yscale("log")
    stitle = "LM Lat (%s) %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)

    figpath = figdir + "/bar.lm-lat.%s.lat%03d-%03d.lon%03d-%03d.png"%(vname,lllat,urlat, lllon, urlon)

plt.suptitle(stitle)
plt.savefig(figpath)
plt.show()


#--- draw (line PDF)-------
fig = plt.figure(figsize=(6,6))
for i,density in enumerate([True, False]):
    ax  = fig.add_subplot(2,1,i+1)
    if density==True:
        a2his = a2rfreq_his
        a2nat = a2rfreq_nat

    else:
        a2his = a2count_his
        a2nat = a2count_nat

    a1his = a2his.mean(axis=0)
    a1nat = a2nat.mean(axis=0)


    a1hisup = np.percentile(a2his, 95, axis=0)
    a1hislw = np.percentile(a2his, 5, axis=0)
    a1natup = np.percentile(a2nat, 95, axis=0)
    a1natlw = np.percentile(a2nat, 5, axis=0)

    ax.plot(acnt, a1his, "-",linewidth=2,color="red")
    ax.plot(acnt, a1nat, "-",linewidth=2,color="blue")

    ax.plot(acnt, a1hisup, "-", linewidth=0.5, color="red")
    ax.plot(acnt, a1hislw, "-", linewidth=0.5, color="red")
    ax.plot(acnt, a1natup, "-", linewidth=0.5, color="blue")
    ax.plot(acnt, a1natlw, "-", linewidth=0.5, color="blue")

    ax.axvline(avehis, linestyle="-", linewidth=1, color="red")
    ax.axvline(avenat, linestyle="-", linewidth=1, color="blue")

    ax.set_yscale("log")
    stitle = "LM Lat (%s) %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)

    figpath = figdir + "/pdf.lm-lat.%s.lat%03d-%03d.lon%03d-%03d.png"%(vname,lllat,urlat, lllon, urlon)

plt.suptitle(stitle)
plt.savefig(figpath)
plt.show()



# %%
