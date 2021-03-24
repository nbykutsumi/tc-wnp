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
#calcflag = True
calcflag = False

prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(1,5+1)
lens    = range(1,50+1)
#xname   = "dslp"
xname   = "wmaxlw"
yname   = "prc"
radkm   = 500

#lpercent=[0,50,90,99,99.9,100]
#lrp =[0,1,10]
lrp =[0,1,5,10]

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

for scen in ["HPB_NAT","HPB"]:
    if calcflag != True: continue


    d1int = {}
    d1prc = {}
    d1now = {}
    for isub, rname in enumerate(lrname): 
        d1int[rname] = deque([])
        d1prc[rname] = deque([])
        d1now[rname] = deque([])

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


        for isub, rname in enumerate(lrname): 
            a1flag = d2mask[rname].flatten()[a1now]
            a1xvecttmp = ma.masked_where(a1flag !=1, a1xvect).compressed()
            a1yvecttmp = ma.masked_where(a1flag !=1, a1yvect).compressed()
            a1nowtmp   = ma.masked_where(a1flag !=1, a1now).compressed()

            d1int[rname].extend(a1xvecttmp)
            d1prc[rname].extend(a1yvecttmp)
            d1now[rname].extend(a1nowtmp)

    #--- save ---
    for rname in lrname:
        datdir = "/home/utsumi/temp/bams2020/decomp-%dkm-pix"%(radkm)
        util.mk_dir(datdir)
        intpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%(xname,scen,rname,lens[0],lens[-1]) 
        prcpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%(yname,scen,rname,lens[0],lens[-1]) 
        nowpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%("now",scen,rname,lens[0],lens[-1]) 

        d1int[rname] = np.array(d1int[rname])
        d1prc[rname] = np.array(d1prc[rname])
        d1now[rname] = np.array(d1now[rname])
        np.save(intpath, d1int[rname])
        np.save(prcpath, d1prc[rname])
        np.save(nowpath, d1now[rname])
        print(intpath)

#-- histogram --
for rname in lrname:

    if xname=="dslp":
        abnd = np.arange(0,25,0.5)
        acnt= (abnd[:-1] + abnd[1:])*0.5
    elif xname=="wmaxlw":
        abnd = np.arange(12,100,2)
        acnt= (abnd[:-1] + abnd[1:])*0.5


    for scen in lscen:
        datdir = "/home/utsumi/temp/bams2020/decomp-%dkm-pix"%(radkm)
        intpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%(xname,scen,rname,lens[0],lens[-1]) 
        prcpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%(yname,scen,rname,lens[0],lens[-1]) 
        nowpath = datdir + "/%s.%s.%s.%03d-%03d.npy"%("now",scen,rname,lens[0],lens[-1]) 

        if scen=="HPB":
            a1int_his = np.load(intpath)
            a1prc_his = np.load(prcpath)
            a1now_his = np.load(nowpath)
        elif scen=="HPB_NAT":
            a1int_nat = np.load(intpath)
            a1prc_nat = np.load(prcpath)
            a1now_nat = np.load(nowpath)
        else:
            print("check scen",scen)
            sys.exit()



    a1numint_his,_   = np.histogram(a1int_his, bins=abnd)
    a1numint_nat,_   = np.histogram(a1int_nat, bins=abnd)

    a2numext_his = np.full([len(lrp),len(abnd)-1], -9999)
    a2numext_nat = np.full([len(lrp),len(abnd)-1], -9999)
    for irp,rp in enumerate(lrp):
        #-- Load extreme criteria (NAT) ---
        if rp ==0:
            a2th = np.zeros([nyin,nxin],float32)

        else:
            [[lllatext,lllonext],[urlatext,urlonext]] =[[0,100],[50,150]]
            
            yext0 = bisect_left(a1latin, lllatext)
            yext1 = bisect_left(a1latin, urlatext)
            xext0 = bisect_left(a1lonin, lllonext)
            xext1 = bisect_left(a1lonin, urlonext)

            thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/HPB_NAT.1990-2010"
            thpath= thdir + "/prec.rp-%03d.WNP.npy"%(rp)
            a2thTmp = np.load(thpath) *60*60*24    # mm/day 

            a2th = np.full([nyin,nxin],-9999).astype(float32)
            a2th[yext0:yext1+1, xext0:xext1+1] = a2thTmp


        a1th_his = a2th.flatten()[a1now_his]
        a1th_nat = a2th.flatten()[a1now_nat]

        a1mask_his = ma.masked_less(a1prc_his, a1th_his).mask
        a1mask_nat = ma.masked_less(a1prc_nat, a1th_nat).mask

        a1int_his_tmp = ma.masked_where(a1mask_his, a1int_his).compressed()
        a1int_nat_tmp = ma.masked_where(a1mask_nat, a1int_nat).compressed()
        a1prc_his_tmp = ma.masked_where(a1mask_his, a1prc_his).compressed()
        a1prc_nat_tmp = ma.masked_where(a1mask_nat, a1prc_nat).compressed()

        a1numext_his,_,_ = scipy.stats.binned_statistic(a1int_his_tmp, a1prc_his_tmp, statistic="count", bins=abnd)
        a1numext_nat,_,_ = scipy.stats.binned_statistic(a1int_nat_tmp, a1prc_nat_tmp, statistic="count", bins=abnd)

        a2numext_his[irp] = a1numext_his
        a2numext_nat[irp] = a1numext_nat

        print(rp)
        print(rp,a1numint_his.sum(), a1numext_his.sum(),a1numint_nat.sum(),a1numext_nat.sum())

    a2prbext_his = (ma.masked_where(np.broadcast_to(a1numint_his, a2numext_his.shape)==0, a2numext_his) / a1numint_his).filled(0)
    a2prbext_nat = (ma.masked_where(np.broadcast_to(a1numint_nat, a2numext_nat.shape)==0, a2numext_nat) / a1numint_nat).filled(0)

    print(a2prbext_his[0])
    print(a2prbext_nat[0])
    #--
    a2dN = (a1numint_his - a1numint_nat)*a2prbext_nat
    a2dT = (a2prbext_his - a2prbext_nat)*a1numint_nat
    a2dNdT = (a1numint_his - a1numint_nat)*(a2prbext_his - a2prbext_nat)
    a2dnum = a2numext_his - a2numext_nat
    a2sum  = a2dN + a2dT + a2dNdT

    a1numext_nat = a2numext_nat.sum(axis=1)
    a1numext_his = a2numext_his.sum(axis=1)
    a1dN   = a2dN.sum(axis=1)
    a1dT   = a2dT.sum(axis=1)
    a1dNdT = a2dNdT.sum(axis=1)
    a1dnum = a2dnum.sum(axis=1)

    lout  = [["rp"] + lrp, ["num_nat"]+a1numext_nat.tolist(), ["num_his"]+a1numext_his.tolist(), ["dnum"]+a1dnum.tolist(),  ["dN"]+a1dN.tolist(), ["dT"]+a1dT.tolist(), ["dNdT"]+a1dNdT.tolist() ]
    sout  = util.list2csv(lout)

    csvpath = datdir + "/decomp.%s.%s.csv"%(xname,rname)
    f=open(csvpath, "w"); f.write(sout); f.close()
    print(csvpath)

    #-- write a2numprc ---
    hispath = datdir + "/a2numext.%s.%s.%s.csv"%(xname,"HPB",rname)
    natpath = datdir + "/a2numext.%s.%s.%s.csv"%(xname,"HPB_NAT",rname)

    sout = util.array2csv(a2numext_his)
    f=open(hispath, "w"); f.write( sout); f.close()

    sout = util.array2csv(a2numext_nat)
    f=open(natpath, "w"); f.write( sout); f.close()

    print(hispath)


##--- draw -------
#Ncol = 1
#Nrow = 2
#
#Fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharex=True, figsize=(6,6))
#Axs = axs.flatten()
#For isub,rname in enumerate(lrname):
#    print("draw",isub)
#    ax = axs[isub]
#    a2ave_his = ma.masked_invalid(np.array(dave_his[rname]))*60*60*24 # mm/day
#    a2ave_nat = ma.masked_invalid(np.array(dave_nat[rname]))*60*60*24 # mm/day
#
#    a1ave_his = a2ave_his.mean(axis=0)
#    a1ave_nat = a2ave_nat.mean(axis=0)
#    a1std_his = a2ave_his.std(axis=0)
#    a1std_nat = a2ave_nat.std(axis=0)
#
#    a1up_his = a1ave_his + a1std_his
#    a1lw_his = a1ave_his - a1std_his
#    a1up_nat = a1ave_nat + a1std_nat
#    a1lw_nat = a1ave_nat - a1std_nat
#
#    ax.plot(acnt, a1ave_his, "-",color="r", linewidth=2)
#    ax.plot(acnt, a1ave_nat, "-",color="b", linewidth=2)
#
#    ax.fill_between(acnt, a1lw_his, a1up_his, color="tomato", alpha=0.5)
#    ax.fill_between(acnt, a1lw_nat, a1up_nat, color="royalblue", alpha=0.5)
#
#    #ax.plot(acnt, a1up_his, "--",color="r", linewidth=1)
#    #ax.plot(acnt, a1lw_his, "--",color="r", linewidth=1)
#    #ax.plot(acnt, a1up_nat, "--",color="b", linewidth=1)
#    #ax.plot(acnt, a1lw_nat, "--",color="b", linewidth=1)
#
#    stitle = rname
#    ax.set_title(stitle)
#
#    #if xname=="wmaxlw":
#    #    ax.set_xlim([10,30])
#    #    ax.set_ylim([0,30])
#
#Plt.tight_layout(rect=[0, 0, 1, 0.96])  # save space for suptitle
#
#Ssuptitle = "%s vs %s (%s) %04d-%04d ens:%03d-%03d"%(xname, yname, metric, iY, eY, lens[0],lens[-1])
##ssuptitle = ssuptitle + "\nlat:%03d-%03d lon:%03d-%03d"%(lllat,urlat, lllon, urlon)
#Plt.suptitle(ssuptitle)
#
#Figpath = figdir + "/plot.mul.%s.vs.%s.%s.customreg.png"%(xname,yname, metric)
#Plt.savefig(figpath)
#
#Figpath = figdir + "/plot.mul.%s.vs.%s.%s.customreg.pdf"%(xname,yname, metric)
#Plt.savefig(figpath)
#Plt.show()
#Print(figpath)
#
#
## %%
#
## %%
