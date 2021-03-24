# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
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
import scipy.stats
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 1990,1990
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
lYM = util.ret_lYM([iY,1],[eY,12])

#-----------------
figflag = True

prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,10+1)
lens    = range(1,50+1)
#lens    = [1]
vname   = "dslp"
#vname   = "wmaxlw"
#vname   = "lat"

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

#--- custom region masks ----
lrname = ["MD","ST"]
d2mask = {}
for rname in lrname:
    maskdir = "/home/utsumi/temp/bams2020/mask"
    maskpath= maskdir + "/mask.%s.npy"%(rname)
    d2mask[rname] = np.load(maskpath)

##-- Land sea mask --
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#a2land = a2land[y0:y1+1, x0:x1+1]

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
dvect_his = {}
dvect_nat = {}

for scen in ["HPB","HPB_NAT"]:
    if vname=="dslp":
        abnd = np.arange(0,25,0.5)
        acnt= (abnd[:-1] + abnd[1:])*0.5
    elif vname=="wmaxlw":
        #abnd = np.arange(12,100,2)
        abnd = np.arange(12,80,1)
        acnt= (abnd[:-1] + abnd[1:])*0.5


    a1vectall = deque([])
    a1nowall = deque([])
    for iens, ens in enumerate(lens):
        print("load",scen,ens)
        vecdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, vname, scen, ens)
        nowdir = vectbaseDir + '/%s/%s/%s-%03d'%(slabel, "nowpos", scen, ens)

        a1vect = np.concatenate([np.load(vecdir + "/%s.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(vname,lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1now = np.concatenate([np.load(nowdir + "/nowpos.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear])

        if vname=="dslp":
            a1vect = a1vect*0.01  # hPa

        a1vectall.extend(a1vect)
        a1nowall.extend(a1now)
    a1vectall = np.array(a1vectall)
    a1nowall = np.array(a1nowall)
    a1nowall = a1nowall - 1   # 1,2,3.. (Fortran) --> 0,1,2,...(Python)

    for isub, rname in enumerate(lrname): 
        a1flag = d2mask[rname].flatten()[a1nowall]
        a1vecttmp = ma.masked_where(a1flag !=1, a1vectall).compressed()
        a1count, _ = np.histogram(a1vecttmp, bins=abnd, density=False)
        a1rfreq, _ = np.histogram(a1vecttmp, bins=abnd, density=True)
    
        if scen=="HPB":
            dcount_his[isub] = a1count
            drfreq_his[isub] = a1rfreq
            dave_his[isub] = a1vecttmp.mean()
            dvect_his[isub]= a1vecttmp

        elif scen=="HPB_NAT":
            dcount_nat[isub] = a1count
            drfreq_nat[isub] = a1rfreq
            dave_nat[isub] = a1vecttmp.mean()
            dvect_nat[isub]= a1vecttmp

        else:
            print("check scen",scen)
            sys.exit()

#-- PV ---
for isub,rname in enumerate(lrname):
    tv, pv = scipy.stats.ttest_ind(dvect_his[isub], dvect_nat[isub], equal_var=False, nan_policy="omit")
    avehis = dvect_his[isub].mean()
    avenat = dvect_nat[isub].mean()
    difave = avehis - avenat
    print("")
    print(rname)
    print("Ave_his=%.1f"%(avehis))
    print("Ave_nat=%.1f"%(avenat))
    print("p-value=%.3f"%(pv))


#--- draw -------
for density in ["count"]:
    if figflag != True: continue
    ncol = 1
    nrow = 2

    fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, figsize=(6,6))
    axs = axs.flatten()
    for isub,rname in enumerate(lrname):
        print("draw",isub)
        ax = axs[isub]
        if density=="rfreq":
            a1his = drfreq_his[isub]
            a1nat = drfreq_nat[isub]

        elif density=="count":
            a1his = dcount_his[isub]
            a1nat = dcount_nat[isub]

        else:
            print("check density",density)
            sys.exit()   

        wbar = abnd[1]-abnd[0]
        ax.bar(acnt, a1his, width=wbar, alpha=0.5,color="red")
        ax.bar(acnt, a1nat, width=wbar, alpha=0.5,color="blue")

        ax.axvline(dave_his[isub], linestyle="-", linewidth=1, color="red")
        ax.axvline(dave_nat[isub], linestyle="-", linewidth=1, color="blue")

        ax.set_yscale("log")
        stitle = rname
        ax.set_title(stitle)

        ##-- difference ---
        #ax2 = ax.twinx()
        #a1dif = a1his - a1nat
        #ax2.plot(acnt, a1dif, "-",color="k")
        #ax2.set_yscale("log")

        #-- save figure data -------------
        figdatdir = "/home/utsumi/temp/bams2020/fig/dat/pdf-tc-intensity"
        util.mk_dir(figdatdir)
        csvpath = figdatdir + "/count-tc-intensity.%s.csv"%(rname)
        avepath = figdatdir + "/ave.%s.csv"%(rname)

        sout = "hPa,HPB_NAT,HPB\n"
        for i in range(len(acnt)):
            sout = sout + "%s,%s,%s\n"%(acnt[i],a1nat[i], a1his[i])

        sout = sout.strip()
        f=open(csvpath,"w"); f.write(sout); f.close()

        sout = "HPB_NAT,%s\n"%(dave_nat[isub]) 
        sout = sout + "HPB,%s"%(dave_his[isub]) 
        f=open(avepath,"w"); f.write(sout); f.close()
        print(csvpath) 

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # save space for suptitle

    ssuptitle = "%s %04d-%04d ens:%03d-%03d"%(vname, iY, eY, lens[0],lens[-1])
    #ssuptitle = ssuptitle + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
    plt.suptitle(ssuptitle)

    figpath = figdir + "/hist-bar.mul.%s.%scustomreg.png"%(vname,density)
    plt.savefig(figpath)

    figpath = figdir + "/hist-bar.mul.%s.%s.customreg.pdf"%(vname,density)
    plt.savefig(figpath)

    plt.show()



# %%
