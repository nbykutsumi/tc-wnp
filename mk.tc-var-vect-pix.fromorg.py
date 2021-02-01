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
import calendar
from detect_fsub import *
#--------------------------------------
calcobj= True
#calcobj= False
#figflag = True
figflag = False

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
#lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,1+1)
#lens    = range(36,50+1)
lens    = range(1,50+1)
#lens    = range(16,50+1)
#lens    = range(2,15+1)
res     = "320x640"
vname = "prc0500"  # calc from map

detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))  # test
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))

hostname=socket.gethostname()
if hostname=="shui":
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=="well":
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'
else:
    print("check hostname",hostname)
    sys.exit()

outbaseDir = '/home/utsumi/temp/bams2020/vect-tc-var'
util.mk_dir(outbaseDir)
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

#************************************
# d4PDF (Objective detection)
#************************************
#lthsst  = [27,28]
#lthsst  = [27,27.5,28]
lthsst  = [27]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
#lthwcore= [-1,0,1]
lthdura = [36]
#lthwind = [10,13,15]
lthwind = [12]
#lthwdif = [-9999]
lthwdif = [-9999]



lkey = [[thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif]
        for thsst   in lthsst
        for exrvort in lexrvort
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thdura in lthdura
        for thwind in lthwind
        for thwdif in lthwdif
        ]


for (thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif) in lkey:
    const  = ConstCyclone.Const(prj=prj, model=model)
    const['Lat'] = d4PDF.Lat()
    const['Lon'] = d4PDF.Lon()

    const['thsst']   = thsst + 273.15   # K
    const['exrvort'] = exrvort
    const['tcrvort'] = tcrvort 
    const['thwcore'] = thwcore
    const['thdura']  = thdura
    const['thwind']  = thwind
    const['thwdif']  = thwdif

    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    for scen in lscen:
        for ens in lens:

            for Year in lYear:
                if calcobj is not True: continue

                ax = deque()
                ay = deque()
                awmax  = deque()
                adslp  = deque()
                aprc   = deque()

                for Mon in lMon:

                    eday   = calendar.monthrange(Year,Mon)[1]
                    iDTime = datetime(Year,Mon,1,0)
                    eDTime = datetime(Year,Mon,eday,18)
                    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=6))

                    ltclonlat = []
                    print(scen,ens,Year,Mon)
                    #-------------------
                    wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)

                    cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)

                    #-- dslp ---
                    _, dtcxy_slp = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname="slp")
                    _, dtcxy_slpave = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname="slp_mean_adj")

                    #-- wmaxlw --
                    _, dtcxy_wmax = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname="wmaxlw")

                    #-- precip --
                    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")
                    d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)                      
                    a3prec_day = d4sfc.load_6hr_mon("PRECIPI", scen, ens, Year, Mon).reshape(-1,4,nyin,nxin).mean(axis=1)

                    ldtime = dtcxy_wmax.keys()

                    for dtime in ldtime:
                        doy = (dtime - datetime(Year,Mon,1)).days  # 0, 1, 2, ...
                        lxywmax = dtcxy_wmax[dtime]
                        if len(lxywmax)==0: continue

                        lwmax   = np.array(lxywmax)[:,2]
                        lslp    = np.array(dtcxy_slp[dtime])[:,2]
                        lslpave = np.array(dtcxy_slpave[dtime])[:,2]
                        if len(lwmax)==0: continue

                        for i,(x,y,_) in enumerate(lxywmax):

                            #print(dtime,x,y)
                            if (x<x0)or(x>x1)or(y<y0)or(y>y1): continue

                            a2mask = detect_fsub.mk_a2mask_with_table(a2table.T, [x], [y], nxin, nyin).T

                            aprcTmp = ma.masked_where(a2mask==0, a3prec_day[doy]).compressed()
                            nrec  = len(aprcTmp)
                            axTmp = np.full(nrec, x)
                            ayTmp = np.full(nrec, y)
                            awmaxTmp = np.full(nrec, lwmax[i])
                            adslpTmp = np.full(nrec, lslp[i]-lslpave[i])

                            ax.extend(axTmp)
                            ay.extend(ayTmp)
                            awmax.extend(awmaxTmp)
                            adslp.extend(adslpTmp)
                            aprc.extend(aprcTmp)

                if len(ax)==0: continue
                ax = np.array(ax).astype("int16")
                ay = np.array(ay).astype("int16")
                awmax = (np.array(awmax)*100).astype("int16")  # scaled by 100
                adslp = np.array(adslp).astype("int16")
                aprc  = (np.array(aprc)*60*60*24).astype("int16")   # mm/day

                #-- save annual data -----
                outDir = outbaseDir + '/%s/tc-vect-500km-pix/%s-%03d'%(slabel, scen, ens)
                util.mk_dir(outDir)
                xpath = outDir + "/x.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)
                ypath = outDir + "/y.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)
                wmaxpath = outDir + "/wmaxlw.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)
                dslppath = outDir + "/dslp.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)
                prcpath = outDir + "/prc.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(lllat,urlat,lllon,urlon,Year)

                np.save(xpath, ax)
                np.save(ypath, ay)
                np.save(wmaxpath, awmax)
                np.save(dslppath, adslp)
                np.save(prcpath, aprc)

                print(prcpath)
#####*** Draw PDF *************
#for scen in ["HPB","HPB_NAT"]:
#    if figflag !=True: continue
#
#    if vname=="dslp":
#        cmbnd = np.arange(0,25,0.5)
#        cmcnt= (cmbnd[:-1] + cmbnd[1:])*0.5
#    elif vname=="wmaxlw":
#        cmbnd = np.arange(12,100,2)
#        cmcnt= (cmbnd[:-1] + cmbnd[1:])*0.5
#
#
#    a2count = np.zeros([len(lens), len(cmcnt)], int32)
#    a2rfreq = np.zeros([len(lens), len(cmcnt)], float64)
#
#    a1vectall = deque([])
#    for iens, ens in enumerate(lens):
#        outDir = outbaseDir + '/%s/%s/%s-%03d'%(slabel, vname, scen, ens)
#        util.mk_dir(outDir)
#
#        a1vect = np.concatenate([np.load(outDir + "/%s.lat%03d-%03d.lon%03d-%03d.%04d.npy"%(vname,lllat,urlat,lllon,urlon,Year)) for Year in lYear])
#
#        if vname=="dslp":
#            a1vect = a1vect*0.01  # hPa
#        
#        a1vectall.extend(a1vect)
#
#        a1count, _ = np.histogram(a1vect, bins=cmbnd, density=False)
#        a1rfreq, _ = np.histogram(a1vect, bins=cmbnd, density=True)
#
#        a2count[iens] = a1count 
#        a2rfreq[iens] = a1rfreq
#    
#    if scen=="HPB":
#        a2count_his = a2count
#        a2rfreq_his = a2rfreq
#        avehis = np.array(a1vectall).mean()
#
#    elif scen=="HPB_NAT":
#        a2count_nat = a2count
#        a2rfreq_nat = a2rfreq
#        avenat = np.array(a1vectall).mean()
#
#    else:
#        print("check scen",scen)
#        sys.exit()
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
#    figpath = figdir + "/pdf.%s.lat%03d-%03d.lon%03d-%03d.png"%(vname,lllat,urlat, lllon, urlon)
#
#plt.suptitle(stitle)
#plt.savefig(figpath)
#plt.show()
#


# %%
