# %%
import matplotlib
%matplotlib inline

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from numpy import *
import numpy as np
import d4PDF
import os, sys
import myfunc.util as util
import socket
from bisect import bisect_left

hostname=socket.gethostname()
if hostname=="shui":
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=="well":
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'
else:
    print("check hostname",hostname)
    sys.exit()

outdir = "/home/utsumi/temp/bams2020/env-ts"
util.mk_dir(outdir)

#calcflag = True
calcflag = False
figflag  = True

region = "WNP"
dbbox = {"WNP":[[0,100],[50,180]]}
[[lllat,lllon],[urlat,urlon]] = dbbox[region]

a1lat = d4PDF.Lat()
a1lon = d4PDF.Lon()
ny = len(a1lat)
nx = len(a1lon)

x0 = bisect_left(a1lon, lllon)
x1 = bisect_left(a1lon, urlon)
y0 = bisect_left(a1lat, lllat)
y1 = bisect_left(a1lat, urlat)

a1latreg = a1lat[y0:y1+1]
a1lonreg = a1lon[x0:x1+1]
nyreg    = len(a1latreg)
nxreg    = len(a1lonreg)

iYear,eYear = 1990,2010
lYear = range(iYear,eYear+1)
lMon  = range(1,12+1)
#lens  = range(1,50+1)
lens  = range(1,7+1)
lscen = ["HPB","HPB_NAT"]
#lscen = ["HPB"]
#lscen = ["HPB_NAT"]

d4 = d4PDF.avr_mon_320x640(vtype="sfc", dbbaseDir=d4pdfdir)

#--- load grid area data ---
areapath = "/home/utsumi/bin/tc-wnp/gridarea.320x640.npy"
a2area = np.load(areapath)[y0:y1+1, x0:x1+1]
a2wgt  = a2area / a2area.sum()
#---------------------------
#lvname = ["TA","PWATER"]
#lvname = ["TA"]
lvname = ["PWATER"]
for vname in lvname:
    if calcflag != True: continue
    for scen in lscen:
        for ens in lens:
            a1v = []
            for Year in lYear:
                a3var = np.array([d4.load_ave_mon(vname, scen, ens, Year, Mon)[y0:y1+1,x0:x1+1] for Mon in lMon])
                v = (a2wgt * a3var.mean(axis=0)).sum()
                a1v.append(v)
                print(scen,ens,Year,v)

            a1v = np.array(a1v).astype("float64") 
            #-- save ---
            tcpath = outdir + "/ts.%s.%s.%03d.%04d-%04d.%s.bn"%(vname,scen,ens,iYear,eYear,region)
            a1v.tofile(tcpath)
            print(tcpath) 

#---- Draw time series -------
for vname in lvname:
    if figflag != True: continue

    #-- load ---
    a2his = np.array([np.fromfile( outdir + "/ts.%s.%s.%03d.%04d-%04d.%s.bn"%(vname,"HPB",ens,iYear,eYear,region), "float64") for ens in lens])

    a2nat = np.array([np.fromfile( outdir + "/ts.%s.%s.%03d.%04d-%04d.%s.bn"%(vname,"HPB_NAT",ens,iYear,eYear,region), "float64") for ens in lens])

    if vname=="TA":
        a2his = a2his - 273.15
        a2nat = a2nat - 273.15


    a1his = a2his.mean(axis=0)
    a1nat = a2nat.mean(axis=0)
    a1x = list(lYear)

    fig = plt.figure(figsize=(6,4))
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])

    ax.plot(a1x, a1his, "-", color="r", linewidth=2)
    ax.plot(a1x, a1nat, "-", color="b", linewidth=2)
    for i in range(len(lens)):
        ax.plot(a1x, a2his[i], "-", color="r", linewidth=0.5)
        ax.plot(a1x, a2nat[i], "-", color="b", linewidth=0.5)

    stitle = "%s %s ens:%03d-%03d"%(vname, region, lens[0],lens[-1])
    plt.title(stitle)
    ax.xaxis.set_major_locator(plt.MultipleLocator(3))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))

    ax.grid()
    ax.grid(which='minor', linestyle='--')
    
# %%
