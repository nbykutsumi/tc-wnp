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
#lens    = range(1,2+1)
lens    = range(1,50+1)
xname   = "dslp"
#xname   = "wmaxlw"
yname  = "prc"
radkm  = 500

#metric = "mean"
#metric = "p99"
metric = "p90"
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

#----------------------------
def percentile99(x):
    return np.nanpercentile(x, 99)

def percentile90(x):
    return np.nanpercentile(x, 90)



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

dave_his = {}
dave_nat = {}
#dstd_his = {}
#dstd_nat = {}
for isub, rname in enumerate(lrname): 
    dave_his[rname] = [] 
    dave_nat[rname] = [] 
    #dstd_his[rname] = [] 
    #dstd_nat[rname] = [] 

for scen in ["HPB","HPB_NAT"]:
    if xname=="dslp":
        abnd = np.arange(0,25,0.5)
        acnt= (abnd[:-1] + abnd[1:])*0.5
    elif xname=="wmaxlw":
        abnd = np.arange(12,100,2)
        acnt= (abnd[:-1] + abnd[1:])*0.5


    for iens, ens in enumerate(lens):
        print("load",scen,ens)
        vectdir = vectbaseDir   + '/%s/tc-vect-%03dkm-pix/%s-%03d'%(slabel, radkm, scen, ens)

        a1locx  = np.concatenate([np.load(vectdir + "/x.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear]).astype("int32")
        a1locy  = np.concatenate([np.load(vectdir + "/y.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)) for Year in lYear]).astype("int32")

        a1xvect = np.concatenate([np.load(vectdir + "/%s.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(xname,lllat,urlat,lllon,urlon,Year)) for Year in lYear])
        a1yvect = np.concatenate([np.load(vectdir + "/%s.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(yname,lllat,urlat,lllon,urlon,Year)) for Year in lYear])

        a1now = a1locx + nxin*a1locy  # 0,1,2, .... nx*ny-1 (Python)

        if xname=="dslp":
            a1xvect = a1xvect.astype("float32")*(-0.01)  # hPa, calced as center - env, for some reason..
        elif xname=="wmaxlw":
            a1xvect = a1xvect.astype("float32")*0.01  # saved with scaled by 100

        a1yvect = a1yvect.astype("float32")   # mm/day

        for isub, rname in enumerate(lrname):
            a1flag = d2mask[rname].flatten()[a1now]
            a1xvecttmp = ma.masked_where(a1flag !=1, a1xvect).compressed()
            a1yvecttmp = ma.masked_where(a1flag !=1, a1yvect).compressed()


            if metric=="mean":
                a1ave,_,_ = scipy.stats.binned_statistic(a1xvecttmp, a1yvecttmp, statistic="mean", bins=abnd)
            elif metric=="p90":
                a1ave,_,_ = scipy.stats.binned_statistic(a1xvecttmp, a1yvecttmp, statistic=percentile90, bins=abnd)
            elif metric=="p99":
                a1ave,_,_ = scipy.stats.binned_statistic(a1xvecttmp, a1yvecttmp, statistic=percentile99, bins=abnd)

            else:
                print("check metric", metric)
            #a1std,_,_ = scipy.stats.binned_statistic(a1xvecttmp, a1yvecttmp, statistic="mean", bins=abnd)

            if scen=="HPB":
                dave_his[rname].append(a1ave)
                #dstd_his[rname].append(a1std)

            elif scen=="HPB_NAT":
                dave_nat[rname].append(a1ave)
                #dstd_nat[rname].append(a1std)

            else:
                print("check scen",scen)
                sys.exit()

#--- draw -------
ncol = 1
nrow = 2

fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, figsize=(6,6))
axs = axs.flatten()
for isub,rname in enumerate(lrname):
    print("draw",isub)
    ax = axs[isub]
    a2ave_his = ma.masked_invalid(np.array(dave_his[rname])) # mm/day
    a2ave_nat = ma.masked_invalid(np.array(dave_nat[rname])) # mm/day

    a1ave_his = a2ave_his.mean(axis=0)
    a1ave_nat = a2ave_nat.mean(axis=0)
    a1std_his = a2ave_his.std(axis=0)
    a1std_nat = a2ave_nat.std(axis=0)

    a1up_his = a1ave_his + a1std_his
    a1lw_his = a1ave_his - a1std_his
    a1up_nat = a1ave_nat + a1std_nat
    a1lw_nat = a1ave_nat - a1std_nat

    ax.plot(acnt, a1ave_his, "-",color="r", linewidth=2)
    ax.plot(acnt, a1ave_nat, "-",color="b", linewidth=2)

    ax.fill_between(acnt, a1lw_his, a1up_his, color="tomato", alpha=0.5)
    ax.fill_between(acnt, a1lw_nat, a1up_nat, color="royalblue", alpha=0.5)

    #ax.plot(acnt, a1up_his, "--",color="r", linewidth=1)
    #ax.plot(acnt, a1lw_his, "--",color="r", linewidth=1)
    #ax.plot(acnt, a1up_nat, "--",color="b", linewidth=1)
    #ax.plot(acnt, a1lw_nat, "--",color="b", linewidth=1)

    stitle = rname
    ax.set_title(stitle)

    #if xname=="wmaxlw":
    #    ax.set_xlim([10,30])
    #    ax.set_ylim([0,30])

    #-- save figure data---------
    figdatdir = "/home/utsumi/temp/bams2020/fig/dat/precip-for-%s"%(xname)
    csvpath = figdatdir + "/prec.%s.%s.%s.csv"%(metric, scen, rname)

    sout = "TC-intensity,HPB_mean,HPB_low,HPB_up,HPB_NAT_mean,HPB_NAT_low,HPB_NAT_up\n"
    for i in range(len(acnt)):
        sout = sout + "%s,%s,%s,%s,%s,%s,%s\n"%(acnt[i],a1ave_his[i],a1lw_his[i],a1up_his[i],a1ave_nat[i],a1lw_nat[i],a1up_nat[i])
    sout=sout.strip()

    util.mk_dir(figdatdir)
    f=open(csvpath, "w"); f.write(sout); f.close()
    print(csvpath)
    #----------------------------

plt.tight_layout(rect=[0, 0, 1, 0.96])  # save space for suptitle

ssuptitle = "%s vs %s pix (%s) %04d-%04d ens:%03d-%03d"%(xname, yname, metric, iY, eY, lens[0],lens[-1])
#ssuptitle = ssuptitle + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
plt.suptitle(ssuptitle)

figpath = figdir + "/plot.mul.%s.vs.%s.pix.%s.customreg.png"%(xname,yname, metric)
plt.savefig(figpath)

figpath = figdir + "/plot.mul.%s.vs.%s.pix.%s.customreg.pdf"%(xname,yname, metric)
plt.savefig(figpath)
plt.show()
print(figpath)


# %%

# %%
