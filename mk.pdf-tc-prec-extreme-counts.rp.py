# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import scipy.stats
#----------------------------------
import sys, os, pickle, socket
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
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
nboot = 5000
#landsea = "land"
landsea = "all"
rp     = 1 # 1, 5, 10, 20
#rp     = 10 # 1, 5, 10, 20
#radkm  = 200 # km
radkm  = 500 # km
figdir  = '/home/utsumi/temp/bams2020/fig/tc-prec'
util.mk_dir(figdir)

hostname=socket.gethostname()
if hostname =='shui':
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'


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

#-- Land sea mask --
a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
a2land = a2land[y0:y1+1, x0:x1+1]
a2seamask = ma.masked_equal(a2land, 0).mask
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
            a2numTmp  = a3numTmp.mean(axis=0) * rp
            a3num[iens] = a2numTmp

    #-- bootstrap --
    if bootstrap==True:
        random.seed(time.time())
        lseq = list(range(len(lens)*len(lYear)))
        for i in range(nboot):
            #li = random.choices(lseq, k=len(lens)*len(lYear))  # resample 1050-yr from 1050-yr data
            li = random.choices(lseq, k=len(lYear))  # resample 20-yr from 1050-yr data
            a3numTmp = a3pool[li]
            a2numTmp = a3numTmp.mean(axis=0) * rp
            a3num[i] = a2numTmp

    if scen=="HPB":
        a3dat_his = a3num
    elif scen=="HPB_NAT":
        a3dat_nat = a3num

    if bootstrap==True:
        if scen=="HPB":
            a3pool_his = a3pool
        elif scen=="HPB_NAT":
            a3pool_nat = a3pool

#***********************
# Figure: ensemble mean
#***********************
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6,6))
axs = axs.flatten()

for isub,subregion in enumerate(lrname):
    ax = axs[isub]

    if landsea=="land":
        a2mask = d2mask[subregion] + a2seamask
    else:
        a2mask = d2mask[subregion]

    a1his = ma.masked_where(np.broadcast_to(a2mask, a3dat_his.shape), a3dat_his).mean(axis=(1,2))
    a1nat = ma.masked_where(np.broadcast_to(a2mask, a3dat_nat.shape), a3dat_nat).mean(axis=(1,2))

    if rp==1:
        if bootstrap==True:
            if landsea=="all":
                if subregion=="MD":
                    a1bnd = np.arange(0.18,0.7,0.01)
                elif subregion=="ST":
                    a1bnd = np.arange(0.18,0.7,0.01)
            elif landsea=="land":
                if subregion=="MD":
                    a1bnd = np.arange(0.08,0.58,0.01)
                elif subregion=="ST":
                    a1bnd = np.arange(0.08,0.58,0.01)
        else:
            a1bnd = np.arange(0,0.68,0.02)
    elif rp==10:
        if bootstrap==True:
            if landsea=="all":
                if subregion=="MD":
                    a1bnd = np.arange(0.1,1.6,0.03)
                elif subregion=="ST":
                    a1bnd = np.arange(0.1,1.6,0.03)

            elif landsea=="land":
                if subregion=="MD":
                    a1bnd = np.arange(0.08,1.3,0.03)
                elif subregion=="ST":
                    a1bnd = np.arange(0.08,1.3,0.03)


        else:
            a1bnd = np.arange(0,0.68,0.005)


    #--------------------------
    a1cnt = 0.5*(a1bnd[:-1]+a1bnd[1:]) 
    #a1freq_his, _ = np.histogram(a1his, bins=a1bnd, density=True)
    #a1freq_nat, _ = np.histogram(a1nat, bins=a1bnd, density=True)

    a1freq_his = np.histogram(a1his, bins=a1bnd, density=False)[0] / len(a1his)
    a1freq_nat = np.histogram(a1nat, bins=a1bnd, density=False)[0] / len(a1nat)


    ax.plot(a1cnt, a1freq_his, linewidth=3, color="r")
    ax.plot(a1cnt, a1freq_nat, linewidth=3, color="blue")

    avehis =a1his.mean()
    avenat =a1nat.mean()
    ax.axvline(avehis, linestyle="--", linewidth=1.5, alpha=0.9, color="r")
    ax.axvline(avenat, linestyle="--", linewidth=1.5, alpha=0.9, color="blue")

    #ax.set_yscale("log")
    #-- statistical test ---
    tv, pv = scipy.stats.ttest_ind(a1his, a1nat, equal_var=False, nan_policy="omit")
    print(subregion,tv,pv)

    #--- FAR -----
    if bootstrap==True:
        #nhis = ma.masked_where(np.broadcast_to(a2mask, a3pool_his.shape), a3pool_his).sum()
        #nnat = ma.masked_where(np.broadcast_to(a2mask, a3pool_nat.shape), a3pool_nat).sum()
        nhis = ma.masked_where(np.broadcast_to(a2mask, a3dat_his.shape), a3dat_his).sum()
        nnat = ma.masked_where(np.broadcast_to(a2mask, a3dat_nat.shape), a3dat_nat).sum()

    else:
        nhis = ma.masked_where(np.broadcast_to(a2mask, a3dat_his.shape), a3dat_his).sum()
        nnat = ma.masked_where(np.broadcast_to(a2mask, a3dat_nat.shape), a3dat_nat).sum()
    far = 1 - float(nnat)/float(nhis)

    #-- Title ---    
    stitle = "%s (%s) d4PDF rp:%d-year pv=%.2f FAR:%.2f [times/%d-year]"%(subregion, landsea, rp, pv, far, rp)
    ax.set_title(stitle)

    #--- save data ------------
    figdatdir  = '/home/utsumi/temp/bams2020/fig/dat/pdf-tc-extreme-precip'
    util.mk_dir(figdatdir)
    hispath = figdatdir + "/count.%s.rp-%02d-yr.%s.%s.npy"%("HPB",rp,subregion,landsea)
    natpath = figdatdir + "/count.%s.rp-%02d-yr.%s.%s.npy"%("HPB_NAT",rp,subregion, landsea)
    csvpath = figdatdir + "/histogram.rp-%02d-yr.%s.%s.csv"%(rp,subregion, landsea)

    sout = "Num,HPB_NAT,HPB\n"
    for i in range(len(a1cnt)):
        sout = sout+"%s,%s,%s\n"%(a1cnt[i], a1freq_nat[i], a1freq_his[i])

    sout=sout.strip()
    np.save(hispath, a1his.filled(miss))
    np.save(natpath, a1nat.filled(miss))
    f=open(csvpath,"w"); f.write(sout); f.close()
    print(hispath)


plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
#print(a1his)
#print("")
#print(a1nat)


    

# %%
