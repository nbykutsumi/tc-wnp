# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import myfunc.util as util
from bisect import bisect_left
from detect_fsub import *
import socket
import IBTrACS
import APHRODITE
#--------------------------------------
#calcflag = True
calcflag = False
figflag = True
#figflag = False
iY = 1960
#eY = 1980
eY = 2015
#eY = 1995
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
#lMon  = range(7,9+1)

#iYbase,eYbase = 1980,1997
iYbase,eYbase = 1960, 1987
#-----------------
miss_out = -9999.
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

#************************************
radkm = 500  # km
if radkm == 500:
    a2table = np.load("./tab.dydx4mask.APHRO.MA.050deg.nrady-008.0500km.npy")

#lrp = [1,5]  # return period, years
lrp = [10]  # return period, years

#************************************
for rp in lrp:
    if calcflag !=True: continue

    thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/aphro"
    thpath= thdir + "/prec.rp-%03d.%s.%04d-%04d.npy"%(rp, region, iYbase, eYbase)
    a2th  = np.load(thpath) 
    #----------------------------------
    for Year in lYear:
        lDTime = util.ret_lDTime(datetime(Year,1,1),datetime(Year,12,31), timedelta(days=1))

        iDTimeBst = datetime(Year,1,1,0)
        eDTimeBst = datetime(Year,12,31,18)

        #--- best track ---
        dbst   = {}
        bst    = IBTrACS.IBTrACS()
        dlonlat = bst.ret_dlonlat(iDTimeBst,eDTimeBst)

        #--- load precipitation --
        a3prec_day = aph.load_year(Year=Year)[0] # mm/day

        #-------------------------
        a2num = np.zeros([ny,nx],float16)
        for i,DTime in enumerate(lDTime):
            print(DTime)
            _,Mon,Day = DTime.timetuple()[:3]

            lDTimeHours = util.ret_lDTime(datetime(Year,Mon,Day,0), datetime(Year,Mon,Day,18), timedelta(hours=6))

            llonlat = []
            for DTimeHours in lDTimeHours:
                try:
                    llonlat = llonlat + dlonlat[DTimeHours]
                except:
                    pass

            if len(llonlat)==0: continue

            a1lon,a1lat = list(zip(*llonlat))
            a1lon = np.array(a1lon)
            a1lat = np.array(a1lat)

            a1xtmp = ((a1lon - lon0)/dgrid).astype(int32)
            a1ytmp = ((a1lat - lat0)/dgrid).astype(int32)

            a1x = []
            a1y = []
            for y,x in zip(a1ytmp,a1xtmp):
                if (x<0)or(x>=nx)or(y<0)or(y>=ny):continue
                a1x.append(x)
                a1y.append(y)

            if len(a1x)==0: continue

            #-- Make mask ---
            a2mask = detect_fsub.mk_a2mask_with_table(a2table.T, a1x, a1y, nx, ny).T

            #if len(a1x)>5:
            #    plt.imshow(a2mask,origin="lower")
            #    plt.colorbar()
            #    plt.show()
            #    sys.exit()
            #    # %%

            a2mask = ma.masked_less_equal(a2mask,0)
            a2mask = ma.masked_where(a3prec_day[i] <=a2th, a2mask)
            a2mask = ma.masked_where(a3prec_day[i] <0, a2mask)
            a2num  = a2num + a2mask.filled(0).astype("int16")
        #--- Save file ---------
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
        precdir = precbasedir + "/aphro-base-%04d-%04d"%(iYbase,eYbase)
        util.mk_dir(precdir)

        numpath = precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)
        np.save(numpath, a2num.astype("int16"))
        print(numpath)


###*** Draw map *************
for rp in lrp:
    if figflag != True: continue

    lieY = [[1980,1997],[1998,2015]]

    for (iYearFig,eYearFig) in lieY:
        lYearFig = range(iYearFig,eYearFig+1)

        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
        precdir = precbasedir + "/aphro-base-%04d-%04d"%(iYbase,eYbase)

        a2num = np.array([np.load(precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)).astype("int32") for Year in lYearFig]).sum(axis=0)

        a2num = a2num.astype("float32") /len(lYearFig)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)

        #-- Figure --------
        a2fig = a2num

        BBox = [[10,105],[45,150]]
        [[lllat,lllon],[urlat,urlon]] = BBox

        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
        axmap.coastlines()

        gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        xticks = np.arange(-180,180+1,15)
        yticks = np.arange(-90,90+1,15)
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
        axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

        axmap.set_extent([lllon,urlon,lllat,urlat])

        #-- Make new colormap (white at the lower end) --
        upper = matplotlib.cm.jet(np.arange(256))
        #lower = np.ones((int(256/4),4))
        lower = np.ones((int(256/8),4))
        for i in range(3):
            lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
        mycm = np.vstack(( lower, upper ))
        mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

        #-- color boundaries norm ------
        cmbnd = list(np.arange(0,5+0.1,0.5))
        cmlabels = list(np.arange(0,5+1,1))

        #cmbnd = [0,0.5,0.9] + list(np.arange(1.1,5+0.1,0.5))
        #cmlabels = [0,0.5,0.9] + list(np.arange(1.1,5+1,1))


        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
        #--extent and coast lines ------
        axmap.set_extent([lllon,urlon,lllat,urlat])
        axmap.coastlines()

        #--contourf ----------
        X,Y = np.meshgrid(a1lonbnd, a1latbnd)
        im = axmap.pcolormesh(X,Y, a2fig, cmap=mycm, norm=norm)

        #X,Y = np.meshgrid(a1lon, a1lat)
        #im = axmap.scatter(X,Y, c=a2fig, s=2, cmap=mycm, norm=norm)
        #--------------------- 
        cbar = plt.colorbar(im)
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

        #-- title, figure name -----
        figpath = figdir + '/map.prec-tc.count-extreme.aphro.%s.rp-%03d.%04d-%04d.png'%(region,rp,iYearFig,eYearFig)

        stitle = 'Count TC extreme prec. times/%s-yr OBS\n'%(rp) + '%04d-%04d'%(iYearFig,eYearFig)

        axmap.set_title(stitle)


        plt.show()

        plt.savefig(figpath)
        print(figpath)


# %%
