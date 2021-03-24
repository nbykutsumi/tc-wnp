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
#from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar, socket, time
import random
import APHRODITE
#--------------------------------------
#lYear0 = range(1980,1997+1)
#lYear1 = range(1998,2015+1)

lYear0 = range(1960,1987+1)
lYear1 = range(1988,2015+1)

iY0,eY0 = lYear0[0],lYear0[-1]
iY1,eY1 = lYear1[0],lYear1[-1]

iYbase,eYbase = iY0, eY0

hostname = socket.gethostname()
if hostname =='shui':
    dbbaseDir = "/tank/utsumi/data/APHRO/APHRO_V1101"
elif hostname=='well':
    dbbaseDir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101"

region= "MA"
res   = "050"
aph = APHRODITE.v1101(region=region, res=res, dbbaseDir=dbbaseDir)
a1lat = aph.Lat
a1lon = aph.Lon
a1latbnd = aph.LatBnd
a1lonbnd = aph.LonBnd
ny = len(a1lat)  # 140
nx = len(a1lon)  # 180
lon0 = a1lonbnd[0]
lat0 = a1latbnd[0]
dgrid = aph.dLat

figdir  = '/home/utsumi/temp/bams2020/tc-prec'
util.mk_dir(figdir)
radkm = 500 # km
#iYbase,eYbase = 1980,1997
lrp = [1]  # return period, years

#--- custom region masks ----
#lrname = ["NJ","SJ", "KR", "EC","SC", "IC", "PH"]
lrname = ["MD","ST"]
d2mask = {}
for rname in lrname:
    maskdir = "/home/utsumi/temp/bams2020/mask"
    maskpath= maskdir + "/mask.aphro.%s.npy"%(rname)
    d2mask[rname] = ma.masked_equal(np.load(maskpath),0).mask

#**************************
# Load site mask
#--------------------------
thratio = 0.9
datdir = "/home/utsumi/temp/bams2020/aphro"
sitepath = datdir + "/frac.site.%04d-%04d.npy"%(iY0,eY1)
widepath = datdir + "/mask.expand.site.th-%0.1f.%04d-%04d.npy"%(thratio, iY0,eY1)

a2sitemask = np.load(sitepath)
a2sitemask = ma.masked_less(a2sitemask,thratio).mask
##**************************
## Load elevation
##--------------------------
#dbdir = dbbaseDir + "/APHRO_%s/%sdeg"%(region, res)
#elevpath = dbdir + "/gtopo30_050deg.grd"
#a2elev = np.fromfile(elevpath, "float32").reshape(ny,nx)
#
#plt.imshow(ma.masked_less(a2elev,0), origin="lower")
#plt.show()
#sys.exit()
## %%
#**************************
###*** Draw map *************
for rp in lrp:

    ##**************************
    # extreme threshold (for mask)
    ##**************************
    exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/aphro"
    expath= exdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iYbase, eYbase)
    a2ex = np.load(expath)
    #---------------------------

    precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
    precdir = precbasedir + "/aphro-base-%04d-%04d"%(iYbase,eYbase)

    a3num0 = np.array([np.load(precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)).astype("int32") for Year in lYear0])
    a3num1 = np.array([np.load(precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)).astype("int32") for Year in lYear1])


    #-- bootstrap ----
    random.seed(time.time())
    nboot  = 1000
    nyear0 = len(lYear0)
    nyear1 = len(lYear1)
    lseq0 = list(range(nyear0))
    lseq1 = list(range(nyear1))

    a3boot0 = np.zeros([nboot,ny,nx],"float32")*-9999
    a3boot1 = np.zeros([nboot,ny,nx],"float32")*-9999
    for i in range(nboot):
        li0 = random.choices(lseq0, k=nyear0)
        li1 = random.choices(lseq1, k=nyear1)

        a3boot0[i] = a3num0[li0].sum(axis=0)
        a3boot1[i] = a3num1[li1].sum(axis=0)


    a3boot0 = a3boot0.astype("float32")/ len(lYear0) * rp
    a3boot1 = a3boot1.astype("float32")/ len(lYear1) * rp
    #a2num = a2num.astype("float32") /len(lYearFig)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)

    #*************************
    # Figure
    #--------------------------
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6,6))
    axs = axs.flatten()

    for isub, subregion in enumerate(lrname):
        ax = axs[isub]
        a2mask = d2mask[subregion] * a2sitemask
        a1vect0 = ma.masked_where(np.broadcast_to(a2mask, a3boot0.shape), a3boot0).mean(axis=(1,2))
        a1vect1 = ma.masked_where(np.broadcast_to(a2mask, a3boot1.shape), a3boot1).mean(axis=(1,2))

        if rp==1:
            #a1bnd = np.arange(0,0.8,0.04)
            a1bnd = np.arange(0,0.8,0.02)
        elif rp==10:
            a1bnd = np.arange(0,2,0.1)

        a1cnt = 0.5*(a1bnd[:-1]+a1bnd[1:])
        a1freq0, _ = np.histogram(a1vect0, bins=a1bnd, density=True)
        a1freq1, _ = np.histogram(a1vect1, bins=a1bnd, density=True)

        ax.plot(a1cnt, a1freq0, color="blue")
        ax.plot(a1cnt, a1freq1, color="r")

        #-- statistical test ---
        tv, pv = scipy.stats.ttest_ind(a1vect0, a1vect1, equal_var=False, nan_policy="omit")
        print(subregion,tv,pv)
        stitle = "%s OBS rp:%d-year pv=%.2f [times/%d-year]"%(subregion, rp, pv, rp)
        ax.set_title(stitle)


    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
    #print(a1his)

# %%
