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
import GSMaP
import d4PDF
import myfunc.regrid.Regrid as Regrid
#--------------------------------------
#calcflag = True
calcflag = False
figflag = True
#figflag = False

sregrid="UP"

lrp = [1]  # return period, years

iY = 2001
eY = 2010
#eY = 1995
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
#lMon  = range(7,9+1)

iYbase,eYbase = 2001,2010
#-----------------
miss_out = -9999.
hostname = socket.gethostname()
if hostname =='shui':
    dbDir = ""
elif hostname=='well':
    dbDir = "/media/disk2/data/GSMaP"

gs = GSMaP.GSMaP_daily(prj="standard",ver="v6",prdName="daily_G", compressed=True, dbDir=dbDir)

a1latorg = gs.Lat   # 0.1 
a1lonorg = gs.Lon   # 0.1
a1latbndorg = gs.LatBnd
a1lonbndorg = gs.LonBnd

a1latup = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonup = d4PDF.Lon()   # dlon = 0.5625
a1latbndup = d4PDF.LatBnd()
a1lonbndup = d4PDF.LonBnd()

region= "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}   # For data
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

y0 = bisect_left(a1latup, lllat)
y1 = bisect_left(a1latup, urlat)
x0 = bisect_left(a1lonup, lllon)
x1 = bisect_left(a1lonup, urlon)

nyglb = len(a1latup)
nxglb = len(a1lonup)
nyreg = y1-y0+1
nxreg = x1-x0+1

us = Regrid.UpScale()
us(a1latorg, a1lonorg, a1latup, a1lonup, globflag=True)

figdir  = '/home/utsumi/temp/bams2020/tc-prec'
util.mk_dir(figdir)

#************************************
radkm = 500  # km
if radkm == 500:
    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")


#************************************
for rp in lrp:
    if calcflag !=True: continue

    thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/GSMaP.%s.%04d-%04d"%(sregrid,iYbase,eYbase)
    thpath= thdir + "/prec.rp-%03d.%s.npy"%(rp, region)
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

        #-------------------------
        a2num = np.zeros([nyreg,nxreg],float16)
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

            a1x = []
            a1y = []
            for (lon,lat) in llonlat:
                y = bisect_left(a1latbndup, lat)
                x = bisect_left(a1lonbndup, lon)

                if (x<x0)or(x>x1)or(y<y0)or(y>y1):continue
                a1x.append(x)
                a1y.append(y)


            if len(a1x)==0: continue


            #--- load precipitation --
            a2org = gs.load_day(DTime)
            a2prec_day = us.upscale(a2org, pergrid=False, miss_in=-9.9990002e+2, miss_out=-9999)[y0:y1+1,x0:x1+1]

            #-- Make mask ---
            a2mask = detect_fsub.mk_a2mask_with_table(a2table.T, a1x, a1y, nxglb, nyglb).T
            a2mask = a2mask[y0:y1+1,x0:x1+1]

            print(a2mask.max())


            a2mask = ma.masked_less_equal(a2mask,0)
            a2mask = ma.masked_where(a2prec_day <=a2th, a2mask)
            a2mask = ma.masked_where(a2prec_day <0, a2mask)
            a2num  = a2num + a2mask.filled(0).astype("int16")


        #--- Save file ---------
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
        precdir = precbasedir + "/gsmap-%s-base-%04d-%04d"%(sregrid,iYbase,eYbase)
        util.mk_dir(precdir)

        numpath = precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)
        np.save(numpath, a2num.astype("int16"))
        print(numpath)


###*** Draw map *************
for rp in lrp:
    if figflag != True: continue


    precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
    precdir = precbasedir + "/gsmap-%s-base-%04d-%04d"%(sregrid,iYbase,eYbase)

    a2num = np.array([np.load(precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)).astype("int32") for Year in lYear]).sum(axis=0)

    a2num = a2num.astype("float32") /len(lYear)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)

    #-- Figure --------
    a2fig = a2num

    BBox = [[10,105],[45,150]]  # For figure
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

    ##-- Make new colormap (white at the lower end) --
    #upper = matplotlib.cm.jet(np.arange(256))
    ##lower = np.ones((int(256/4),4))
    #lower = np.ones((int(256/8),4))
    #for i in range(3):
    #    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
    #mycm = np.vstack(( lower, upper ))
    #mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

    #-- Make new colormap (white at the lower end) --
    middle = matplotlib.cm.jet(np.arange(256))
    lower = np.ones((int(256/8),4))
    upper = np.ones((int(256/16),4))
    for i in range(3):
        lower[:,i] = np.linspace(1, middle[0,i], lower.shape[0])

    for i in range(3):
        upper[:,i] = np.linspace(middle[-50,i], middle[-20,i], upper.shape[0] )

    mycm = np.vstack(( lower, middle[:-20], upper ))
    mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])
 

    #-- color boundaries norm ------
    cmbnd = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    cmlabels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

    #cmbnd = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    #cmlabels = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]


    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #-- meshgrid ----------
    X,Y = np.meshgrid(a1lonbndup[x0:x1+1], a1latbndup[y0:y1+1])
    im = axmap.pcolormesh(X,Y, a2fig, cmap=mycm, norm=norm)
    cbar = plt.colorbar(im)
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)

    #-- title, figure name -----
    figpath = figdir + '/map.prec-tc.count-extreme.aphro.%s.rp-%03d.%s.%04d-%04d.png'%(region,rp,sregrid,iY,eY)

    stitle = 'Count TC extreme prec. times/%s-yr GSMaP\n'%(rp) + '%04d-%04d %s'%(iY,eY, sregrid)

    axmap.set_title(stitle)


    plt.show()

    plt.savefig(figpath)
    print(figpath)


# %%
