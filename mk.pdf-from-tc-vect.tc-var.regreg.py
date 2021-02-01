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
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 2001,2010
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
#lens    = [1,2,3]
#vname   = "dslp"
vname   = "wmaxlw"
#vname   = "lat"
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
figdir  = '/home/utsumi/temp/bams2020/fig/pdf-tc-var'
util.mk_dir(figdir)
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

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

####*** Draw PDF (No ensemble range)*************

dcount_his = {}
dcount_nat = {}
drfreq_his = {}
drfreq_nat = {}
dave_nat = {}
dave_his = {}
for scen in ["HPB","HPB_NAT"]:
    if vname=="dslp":
        abnd = np.arange(0,25,0.5)
        acnt= (abnd[:-1] + abnd[1:])*0.5
    elif vname=="wmaxlw":
        abnd = np.arange(12,100,2)
        acnt= (abnd[:-1] + abnd[1:])*0.5


    a1vectall = deque([])
    a1latall = deque([])
    a1lonall = deque([])
    for iens, ens in enumerate(lens):
        print("load",scen,ens)
        vecdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, vname, scen, ens)
        latdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "lat", scen, ens)
        londir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "lon", scen, ens)

        a1vect = np.concatenate([np.load(vecdir + "/%s.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(vname,lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1lat = np.concatenate([np.load(latdir + "/lat.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1lon = np.concatenate([np.load(londir + "/lon.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])

        if vname=="dslp":
            a1vect = a1vect*0.01  # hPa

        a1vectall.extend(a1vect)
        a1latall.extend(a1lat)
        a1lonall.extend(a1lon)

    a1vectall = np.array(a1vectall)
    a1latall = np.array(a1latall)
    a1lonall = np.array(a1lonall)

    for isub,subbox in enumerate(lsubbox): 
        [[lllatsub,lllonsub],[urlatsub,urlonsub]] = subbox

        a1latmask = ma.masked_outside(a1latall, lllatsub, urlatsub).mask
        a1lonmask = ma.masked_outside(a1lonall, lllonsub, urlonsub).mask
        a1mask = a1latmask + a1lonmask
        a1vecttmp = ma.masked_where(a1mask, a1vectall).compressed()
        a1count, _ = np.histogram(a1vecttmp, bins=abnd, density=False)
        a1rfreq, _ = np.histogram(a1vecttmp, bins=abnd, density=True)
    
        if scen=="HPB":
            dcount_his[isub] = a1count
            drfreq_his[isub] = a1rfreq
            dave_his[isub] = a1vecttmp.mean()
    
        elif scen=="HPB_NAT":
            dcount_nat[isub] = a1count
            drfreq_nat[isub] = a1rfreq
            dave_nat[isub] = a1vecttmp.mean()

        else:
            print("check scen",scen)
            sys.exit()

#--- draw -------
#for density in ["count","rfreq"]:
for norm in ["org","norm"]:
    for density in ["count"]:
        ncol = 9
        nrow = 8

        if (norm=="norm")&(density=="rfreq"): continue

        if norm=="org":
            fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, figsize=(20,20))
        elif norm=="norm":
            fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, sharey=True, figsize=(20,20))

        axs = axs.flatten()
        for isub,subbox in enumerate(lsubbox):
            print("draw",isub)
            [[lllatsub,lllonsub],[urlatsub,urlonsub]] = subbox
            subregion = "%03d.%03d"%(lllatsub,lllonsub)
            ax = axs[isub]
            #ax  = fig.add_subplot(nline,nrow,isub+1, sharex=True, sharey=True)
            if density=="rfreq":
                a1his = drfreq_his[isub]
                a1nat = drfreq_nat[isub]

            elif density=="count":
                if norm=="norm":
                    natmax = max(dcount_nat[isub].max(), dcount_his[isub].max())
                    a1his = dcount_his[isub] / float(natmax)
                    a1nat = dcount_nat[isub] / float(natmax)
                elif norm=="org":
                    a1his = dcount_his[isub]
                    a1nat = dcount_nat[isub]
                else:
                    print("check norm",norm)
                    sys.exit()

            else:
                print("check density",density)
                sys.exit()   

            wbar = abnd[1]-abnd[0]
            ax.bar(acnt, a1his, width=wbar, alpha=0.5,color="red")
            ax.bar(acnt, a1nat, width=wbar, alpha=0.5,color="blue")

            ax.axvline(dave_his[isub], linestyle="-", linewidth=1, color="red")
            ax.axvline(dave_nat[isub], linestyle="-", linewidth=1, color="blue")

            ax.set_yscale("log")
            stitle = "lat:%03d lon:%03d"%(lllatsub, lllonsub)
            ax.set_title(stitle)

        plt.tight_layout(rect=[0, 0, 1, 0.96])  # save space for suptitle

        ssuptitle = "%s %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1]) + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
        plt.suptitle(ssuptitle)

        figpath = figdir + "/hist-bar.mul.%s.%s.%s.png"%(vname,density,norm)
        plt.savefig(figpath)

        figpath = figdir + "/hist-bar.mul.%s.%s.%s.pdf"%(vname,density,norm)
        plt.savefig(figpath)

        plt.show()



# %%
