# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#----------------------------------
import sys, os, pickle
from   numpy import *
import scipy.stats
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
from bisect import bisect_left
from detect_fsub import *
import socket
from collections import deque
#--------------------------------------
#calcflag = True
calcflag = False
figflag = True  
#figflag = False

#meanflag = "ave"
meanflag = "pix"

iY = 1990
eY = 2010

lYear = range(iY,eY+1)
lMon  = range(1,12+1)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
#lens    = list(range(16,20+1))
#lens    = list(range(23,50+1))
#lens    = list(range(36,50+1))
lens    = list(range(1,50+1))
#lens = [1]
region  = "WNP"

dlatsub  = 5
dlonsub  = 5
llatsub0 = np.arange(5,40+0.1,dlatsub)[::-1]
llonsub0 = np.arange(105,145+0.1,dlonsub)

detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'


compbasedir= '/home/utsumi/temp/bams2020/composite'
figdir  = '/home/utsumi/temp/bams2020/fig/map-prec'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640
miss =-9999
miss_int= -9999

#************************************
precbasedir  = "/home/utsumi/temp/bams2020/tc-prec-all"
vectbasedir  = "/home/utsumi/temp/bams2020/tc-prec-all/vect"
bboxpath = precbasedir + "/HPB.001/bbox.%s.npy"%(region)
yxbboxpath = precbasedir + "/HPB.001/yxbbox.%s.npy"%(region)
bbox   = np.load(bboxpath)
yxbbox = np.load(yxbboxpath)

[[lllat,lllon],[urlat,urlon]] = bbox
[[y0,x0],[y1,x1]] = yxbbox
nyreg = y1-y0+1
nxreg = x1-x0+1
a1latreg = a1lat[y0:y1+1]
a1lonreg = a1lon[x0:x1+1]

#--- Make subregion mask for saved daily array ---
d2mask = {}
#for subregion in lsubregion:
lsubbox = [[[latsub0,lonsub0],[latsub0+dlatsub,lonsub0+dlonsub]] for latsub0 in llatsub0 for lonsub0 in llonsub0]

lsubregion = []
for subbox in lsubbox:
    #subbox = dsubbox[subregion]
    [[lllatsub,lllonsub],[urlatsub,urlonsub]] = subbox
    subregion = "%03d.%03d"%(lllatsub,lllonsub)
    lsubregion.append(subregion)

    xsub0 = bisect_left(a1lonreg,lllonsub) 
    xsub1 = bisect_left(a1lonreg,urlonsub) 
    ysub0 = bisect_left(a1latreg,lllatsub) 
    ysub1 = bisect_left(a1latreg,urlatsub) 

    d2mask[subregion] = np.full([nyreg,nxreg], True, dtype=bool)
    d2mask[subregion][ysub0:ysub1+1,xsub0:xsub1+1] = False

    #plt.imshow(d2mask[subregion], origin="lower")
    #plt.colorbar()
    #plt.show()

#-------------------------------------------------

for scen in lscen:
    if calcflag !=True: continue
    for iens, ens in enumerate(lens):
        #print(scen,ens)

        dprec = {}
        for subregion in lsubregion:
            dprec[subregion] = deque([])


        for itmp,Year in enumerate(lYear):
            print(scen,ens,Year)

            precpath = precbasedir + "/%s.%03d/prec-tc.%s.%03d.npy"%(scen,ens,region,Year)

            a3prec = ma.masked_less(np.load(precpath),0)

            for subregion in lsubregion:
                a2mask = d2mask[subregion]

                if meanflag=="ave":
                    a1prec = ma.masked_where(np.broadcast_to(a2mask, a3prec.shape), a3prec).mean(axis=(1,2)).compressed()*60*60*24 # mm/day
                else:
                    a1prec = ma.masked_where(np.broadcast_to(a2mask, a3prec.shape), a3prec).compressed()*60*60*24 # mm/day

                dprec[subregion].extend(a1prec)

        #-- Save ------------
        vectdir = precbasedir + "/vect/%s.d%02d.d%02d.%04d-%04d/%s.%03d"%(meanflag,dlatsub,dlonsub,iY,eY,scen,ens)
        util.mk_dir(vectdir)

        for subregion in lsubregion:
            vectpath = vectdir + "/tc-prec.%s.npy"%(subregion)
            np.save(vectpath, list(dprec[subregion]))
            print(vectpath)


#**********************************
# Draw histogram (regular subregions)
#**********************************
bnd = np.arange(0,1200+1,20)
bndcnt=(bnd[1:] + bnd[:-1])*0.5
print("draw")
for norm in ["org","norm"]:
    if figflag != True: continue
    ncol = 9
    nrow = 8

    if norm=="org":
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(20,20))
    elif norm=="norm":
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol, sharey=True, figsize=(20,20))

    axs = axs.flatten()
    for isub,subregion in enumerate(lsubregion):
        ax = axs[isub]

        for scen in lscen:
            print("plot",subregion,scen)
            a1vect = []
            for ens in lens:
                vectdir = precbasedir + "/vect/%s.d%02d.d%02d.%04d-%04d/%s.%03d"%(meanflag,dlatsub,dlonsub,iY,eY,scen,ens)

                vectpath = vectdir + "/tc-prec.%s.npy"%(subregion)
                a1vect.extend(np.load(vectpath))

            a1freq, _ = np.histogram(a1vect, bins=bnd, density=False)
            if scen=="HPB":
                a1his = a1freq
            elif scen=="HPB_NAT":
                a1nat = a1freq
            else:
                print("check scen",scen)
                sys.exit()

        slllat,slllon = subregion.split(".")
        stitle = "lat:%s lon:%s"%(slllat,slllon)
        ax.set_title(stitle)

        wbar = bnd[1] - bnd[0]
        ax.bar(bndcnt, a1his, width=wbar, alpha=0.5, color="red", log=True)
        ax.bar(bndcnt, a1nat, width=wbar, alpha=0.5, color="blue", log=True)
        #ax.plot(bndcnt, a1his, color="red")
        #ax.plot(bndcnt, a1nat, color="blue")
        #ax.set_yscale("log")

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # save space for suptitle

    ssuptitle = "precip rate [mm/day] %04d-%04d ens:%03d-%03d"%(iY, eY, lens[0],lens[-1])
    plt.suptitle(ssuptitle)
    

    figdir = "/home/utsumi/temp/bams2020/fig/pdf-tc-var"
    util.mk_dir(figdir)

    figpath= figdir + "/hist.mul.tc-prec.%s.png"%(norm)
    plt.savefig(figpath)

    figpath= figdir + "/hist.mul.tc-prec.%s.pdf"%(norm)
    plt.savefig(figpath)

    print(figpath)
    plt.show()

