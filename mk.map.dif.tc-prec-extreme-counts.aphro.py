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
lYear0 = range(1980,1997+1)
lYear1 = range(1998,2015+1)

#lYear0 = range(1960,1987+1)
#lYear1 = range(1988,2015+1)

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
lrp = [1,5]  # return period, years


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
    n = 500
    nyear0 = len(lYear0)
    nyear1 = len(lYear1)
    lseq0 = list(range(nyear0))
    lseq1 = list(range(nyear1))

    a3boot0 = np.zeros([n,ny,nx],"int32")*-9999
    a3boot1 = np.zeros([n,ny,nx],"int32")*-9999
    for i in range(n):
        li0 = random.choices(lseq0, k=nyear0)
        li1 = random.choices(lseq1, k=nyear1)

        a3boot0[i] = a3num0[li0].sum(axis=0)
        a3boot1[i] = a3num1[li1].sum(axis=0)


    a3boot0 = a3boot0.astype("float32")
    a3boot1 = a3boot1.astype("float32")
    #a2num = a2num.astype("float32") /len(lYearFig)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)
        
    plt.figure()
    plt.imshow(ma.masked_equal(a3boot0,0).mean(axis=0), origin="lower")
    plt.colorbar()
    plt.show()

    plt.figure()
    plt.imshow(ma.masked_equal(a3boot1,0).mean(axis=0), origin="lower")
    plt.colorbar()
    plt.show()

    #ytmp,xtmp = 75,120
    ytmp,xtmp = 103,160   # around Kanto
    aboot0 = a3boot0[:,ytmp,xtmp]
    aboot1 = a3boot1[:,ytmp,xtmp]
    atv, apv = scipy.stats.ttest_ind(aboot0, aboot1, axis=0, equal_var=False, nan_policy="omit")

    plt.figure()
    plt.hist(aboot0, 20, histtype="step", facecolor="b")
    plt.hist(aboot1, 20, histtype="step", facecolor="r")
    plt.show()
    #sys.exit()
    ## %%

    a2freq0 = a3num0.sum(axis=0).astype(float32) / nyear0 * rp  # times/rp-year (e.g., rp=10 --> times/10-year)

    a2freq1 = a3num1.sum(axis=0).astype(float32) / nyear1 * rp  # times/rp-year (e.g., rp=10 --> times/10-year)

    a2dif = a2freq1 - a2freq0

    a2tv, a2pv = scipy.stats.ttest_ind(a3boot0, a3boot1, axis=0, equal_var=False, nan_policy="omit")

    a2fig = a2dif
    a2hatch = ma.masked_where(a2pv>0.05, a2pv)
    a2hatch = ma.masked_invalid(a2hatch)

    #--- mask ---------------
    a2fig = ma.masked_where(a2sitemask, a2fig)
    a2hatch = ma.masked_where(a2sitemask, a2hatch)

    ##***********************
    #sys.exit()
    ## %%
    #plt.imshow(a2hatch, origin="lower")
    #plt.colorbar()
    #plt.show()

    a2fig   = ma.masked_where(a2ex<0, a2fig)
    a2hatch = ma.masked_where(a2ex<0, a2hatch)
    #-- title, figure name -----
    stitle = 'diff. TC extreme prec. times/%s-yr\n %04d-%04d to %04d-%04d'%(rp,iY0,eY0,iY1,eY1) 
    figpath = figdir + '/map.dif.prec-tc.count-extreme.obs.rp-%03d.%04d-%04d.to.%04d-%04d.png'%(rp,iY0,eY0,iY1,eY1)
    #---------------------------
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8], projection=ccrs.PlateCarree())
    axmap.set(facecolor="0.8")

    #-- grid lines ---
    BBox = [[10,105],[45,150]]
    [[lllat,lllon],[urlat,urlon]] = BBox

    gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,90, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    #-- set extent and coastlines----
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines(color="gray")

    #-- color boundaries norm --------
    #cmbnd = list(np.arange(-2.5,2.5+0.01,0.25))
    #cmlabels = list(np.arange(-2.5,2.5+0.01,1))
    if rp in [1]:
        #cmbnd = list(np.arange(-0.9, 0.9+0.001, 0.2))
        cmbnd = list(np.arange(-0.9, 0.9+0.001, 0.2))
        cmlabels= list(map(float, ["%.1f"%x for x in cmbnd]))
    if rp in [5]:
        cmbnd = list(np.arange(-1.7, 1.7+0.001, 0.2))
        cmlabels= list(map(float, ["%.1f"%x for x in cmbnd]))

    elif rp in [5]:
        cmbnd = list(np.arange(-0.3, 0.3+0.001, 0.05))
        cmlabels = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3]
    elif rp in [10,20]:
        cmbnd = list(np.arange(-0.5, 0.5+0.001, 0.1))
        cmlabels = [-0.5,-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]

    mycm = 'RdBu_r'
    
    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    
    #-- pcolormesh --
    X,Y = np.meshgrid(a1lonbnd,a1latbnd)
    vmin, vmax = cmbnd[0], cmbnd[-1]
    im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
    print(cmbnd)

    #-- hatch --------------
    #a2hatch = ma.masked_inside(a2hatch, -0.5, 0.5) # should be adjusted considering 
    #plt.pcolor(X, Y, a2hatch, hatch="////", alpha=0.)

    
    X,Y = np.meshgrid(a1lon,a1lat)
    X = ma.masked_where(a2hatch.mask, X)
    Y = ma.masked_where(a2hatch.mask, Y)

    thcolor = {1:0.5, 5:0.9}[rp]
    Xtmp = ma.masked_where(ma.masked_outside(a2fig, -thcolor, thcolor).mask, X)
    Ytmp = ma.masked_where(ma.masked_outside(a2fig, -thcolor, thcolor).mask, Y)
    plt.plot(Xtmp, Ytmp, ".", markersize=1.3, color="0.2")

    Xtmp = ma.masked_where(ma.masked_inside(a2fig, -thcolor, thcolor).mask, X)
    Ytmp = ma.masked_where(ma.masked_inside(a2fig, -thcolor, thcolor).mask, Y)
    plt.plot(Xtmp, Ytmp, ".", markersize=1.3, color="0.8")




    #-- draw colorbar ------
    cax = figmap.add_axes([0.84,0.2,0.02,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)
    #cbar.set_ticks(cbar.ax.get_yticks())
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)

    #-- Tiele -----
    axmap.set_title(stitle)
    #-- Save ------
    plt.savefig(figpath)
    plt.show()
    print(figpath)



# %%
