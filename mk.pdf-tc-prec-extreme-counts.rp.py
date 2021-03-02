# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import scipy.stats
#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import util
import calendar
from bisect import bisect_left
import random, time
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
lYear = range(iY,eY+1)
#-----------------
prj     = "d4PDF"
#lens    = range(1,20+1)
lens    = range(1,50+1)
#lens    = range(1,2+1)
#lens    = range(1,4+1)

#bootstrap = False
bootstrap = True
nboot = 1000

rp     = 1 # 1, 5, 10, 20
#rp     = 10 # 1, 5, 10, 20
#radkm  = 200 # km
radkm  = 500 # km
figdir  = '/home/utsumi/temp/bams2020/fig/tc-prec'
util.mk_dir(figdir)

#----------------------------------
detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
a1latcnt  = d4PDF.Lat()
a1loncnt  = d4PDF.Lon()
a1latbnd  = d4PDF.LatBnd()
a1lonbnd  = d4PDF.LonBnd()
a1lonbnd[0] = 0
a1lonbnd[-1] = 360
miss_int= -9999
ny = len(a1latcnt)
nx = len(a1loncnt) 
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]

region= "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

y0 = bisect_left(a1latcnt, lllat)
y1 = bisect_left(a1latcnt, urlat)
x0 = bisect_left(a1loncnt, lllon)
x1 = bisect_left(a1loncnt, urlon)

nyreg = y1-y0+1    #
nxreg = x1-x0+1

#**************************
miss =-9999
miss_int= -9999
thsst = 27
exrvort = 3*1e-5
tcrvort = 3*1e-5
thwcore  = 0
thwind   = 12
thwdif   = -9999  # not used
thdura   = 36
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

#--- custom region masks ----
#lrname = ["NJ","SJ", "KR", "EC","SC", "IC", "PH"]
lrname = ["MD","ST"]
d2mask = {}
for rname in lrname:
    maskdir = "/home/utsumi/temp/bams2020/mask"
    maskpath= maskdir + "/mask.%s.npy"%(rname)
    d2mask[rname] = ma.masked_equal(np.load(maskpath),0).mask[y0:y1+1,x0:x1+1]

#************************************
# Load data
#************************************
for scen in ["HPB","HPB_NAT"]:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    if bootstrap==True:
        a3pool = np.zeros([len(lens)*len(lYear),nyreg,nxreg], "float32")
        a3num = np.zeros([nboot,nyreg,nxreg], "float32")
    else:
        a3num = np.zeros([len(lens),nyreg,nxreg], "float32")

    for iens,ens in enumerate(lens):
        print(scen,ens)
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)

        a3numTmp = np.array([np.load(precbasedir + "/%s.%03d/num.%s.rp-%03d.%04d.npy"%(scen,ens,region,rp,Year)) for Year in lYear])

        if bootstrap==True:
            a3pool[iens*len(lYear):(iens+1)*len(lYear)] = a3numTmp

        else:
            a2numTmp  = a3numTmp.sum(axis=0) / float(len(lYear)) * rp
            a3num[iens] = a2numTmp

    #-- bootstrap --
    if bootstrap==True:
        random.seed(time.time())
        lseq = list(range(len(lYear)))
        for i in range(nboot):
            li = random.choices(lseq, k=len(lYear))
            a3numTmp = a3pool[li]
            a2numTmp = a3numTmp.sum(axis=0) / float(len(lYear)) * rp
            a3num[i] = a2numTmp

    if scen=="HPB":
        a3dat_his = a3num
    elif scen=="HPB_NAT":
        a3dat_nat = a3num

#***********************
# Figure: ensemble mean
#***********************
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6,6))
axs = axs.flatten()

for isub,subregion in enumerate(lrname):
    ax = axs[isub]

    a2mask = d2mask[subregion]
    
    a1his = ma.masked_where(np.broadcast_to(a2mask, a3dat_his.shape), a3dat_his).mean(axis=(1,2))
    a1nat = ma.masked_where(np.broadcast_to(a2mask, a3dat_nat.shape), a3dat_nat).mean(axis=(1,2))

    if rp==1:
        a1bnd = np.arange(0,0.8,0.04)
    elif rp==10:
        a1bnd = np.arange(0,2,0.1)


    a1cnt = 0.5*(a1bnd[:-1]+a1bnd[1:]) 
    a1freq_his, _ = np.histogram(a1his, bins=a1bnd, density=True)
    a1freq_nat, _ = np.histogram(a1nat, bins=a1bnd, density=True)

    ax.plot(a1cnt, a1freq_his, color="r")
    ax.plot(a1cnt, a1freq_nat, color="blue")
    #-- statistical test ---
    tv, pv = scipy.stats.ttest_ind(a1his, a1nat, equal_var=False, nan_policy="omit")
    print(subregion,tv,pv)
    stitle = "%s d4PDF rp:%d-year pv=%.2f [times/%d-year]"%(subregion, rp, pv, rp)
    ax.set_title(stitle)


plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
#print(a1his)
#print("")
#print(a1nat)


    

# %%
