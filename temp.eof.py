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
from sklearn.decomposition import PCA
from sklearn.decomposition import TruncatedSVD
from sklearn.utils.extmath import randomized_svd
#from scipy.sparse.linalg import svds
#from scipy.linalg import svd as SVD
#--------------------------------------
lYear = range(1980,2010+1)
#-----------------
prj     = "d4PDF"
#lens    = range(1,20+1)
#lens    = range(1,50+1)
lens    = range(1,5+1)
lscen = ["HPB"]
rp     = 1 # 1, 5, 10, 20
#radkm  = 200 # km
radkm  = 500 # km

iYbase = 1990
eYbase = 2010
lensbase = list(range(1,50+1))
scenbase = "HPB_NAT"

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

dbbaseDir = "/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM"
d4sst = d4PDF.avr_mon_320x640(vtype="sfc", dbbaseDir=dbbaseDir)

#************************************
def draw_map(a2dat, stitle, cmbnd):
    a2fig = a2dat
    #---------------------------
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
    axmap.set_facecolor("0.8")
    #-- grid lines ---
    gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
    xticks   = np.arange(-180, 180+1, 15)
    yticks   = np.arange(-90,90, 15)
    gl.xlocator = mticker.FixedLocator(xticks)
    gl.ylocator = mticker.FixedLocator(yticks)
    axmap.set_xticks(xticks, crs=ccrs.PlateCarree())
    axmap.set_yticks(yticks, crs=ccrs.PlateCarree())
    #-- set extent and coastlines----
    #axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.set_extent([105,150,10,45])
    axmap.coastlines(color="k")

    #-- color boundaries norm --------
    cmlabels = cmbnd
    mycm = 'RdBu_r'

    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

    X,Y = np.meshgrid(a1lonbnd[x0:x1+1],a1latbnd[y0:y1+1])
    vmin, vmax = cmbnd[0], cmbnd[-1]
    im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
    print(cmbnd)

    #-- draw colorbar ------
    cax = figmap.add_axes([0.84,0.2,0.02,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)
    #cbar.set_ticks(cbar.ax.get_yticks())
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)

    #-- Tiele -----
    axmap.set_title(stitle)
    #-- Save ------
    return fig

#************************************
# Load IPO index
#------------------------------------
datpath = "./tpi.timeseries.ersstv5.data"
f=open(datpath,"r");lines=f.readlines(); f.close()
lipo = []
for line in lines[1:]:
    line = line.split()
    Year = int(line[0])
    if Year ==-99: break
    if Year in lYear:
        vmean = mean(list(map(float, line[1:])))
        print(vmean)
        lipo.append(vmean)
a1ipo = np.array(lipo)
#************************************
# Load total TC precip (HPB) for mask
#************************************
precbasedir = "/home/utsumi/temp/bams2020/tc-prec-%04dkm"%(radkm)
#avedir = precbasedir + "/ens-ave-%04d-%04d"%(iY,eY)  # this data is created by mk.map.tc-precip.py
avedir = precbasedir + "/ens-ave-%04d-%04d"%(1990,2010)  # this data is created by mk.map.tc-precip.py

a3sum = np.array([np.load(avedir + "/precsum-tc-ave.%s.%03d.npy"%("HPB",ens)) for ens in lens])

a3sum = ma.masked_less(a3sum, 0)

a2tcprec = a3sum.mean(axis=0)[y0:y1+1,x0:x1+1] # unit: mm/year

#************************************
# Load data
#************************************
for scen in lscen:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    a3num = np.zeros([len(lYear),nyreg,nxreg], "float32")
    a3sst = np.zeros([len(lYear),nyreg,nxreg], "float32")
    for i, Year in enumerate(lYear):
        print(Year)
        #precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)

        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm-base-%s-ens-%03d-%03d-yr-%04d-%04d"%(radkm, scenbase, lensbase[0], lensbase[-1], iYbase, eYbase)

        a2numTmp = np.array([np.load(precbasedir + "/%s.%03d/num.%s.rp-%03d.%04d.npy"%(scen,ens,region,rp,Year)) for ens in lens]).mean(axis=0)

        a3num[i] = a2numTmp

        #--- SST ----
        a2sstTmp = np.array([np.load("/home/utsumi/temp/bams2020/map-year/TGEF/HPB.%03d/TGEF.%04d.320x640.npy"%(ens, Year)) for ens in lens]).mean(axis=0)[y0:y1+1, x0:x1+1]
        a3sst[i] = a2sstTmp

    #*********************************
    # SVD
    #---------------------------------
    a2num = a3num.reshape(len(lYear),-1)
    a2sst = a3sst.reshape(len(lYear),-1)
    
    a2num = a2num - a2num.mean(axis=0)
    a2sst = a2sst - a2sst.mean(axis=0)

    a2cov = np.matmul(a2num.T, a2sst)

    U,S,VT = randomized_svd(a2cov, n_components=3)
    V = VT.T
    print("U1",U)
    print("V1",V)

    #----
    tsvd = TruncatedSVD(n_components=3)
    tsvd.fit(a2cov)
    V2 = tsvd.components_.T
    print("trancated",tsvd.explained_variance_ratio_)
    print("V2",V2)
    print("V2.shape",V2.shape)
    #----

    TU = np.matmul(a2num, U)
    TV = np.matmul(a2sst, V)
    print(TU.shape,TV.shape)
    for i in range(TU.shape[1]):
        fig = plt.figure()
        ax  = fig.add_subplot()
        ax.plot(lYear,TU[:,i], "-",label="TU%d"%(i))
        ax.plot(lYear,TV[:,i], "--",label="TV%d"%(i))

        ax2 = ax.twinx()
        #ax2.plot(lYear,a3sst.mean(axis=(1,2)), color="r")
        if i==0:
            ax2.plot(lYear,-a1ipo, "--",color="r")
        else:
            ax2.plot(lYear,a1ipo, "-",color="r")

        ax.legend()
        plt.show()
    #%%
    ##*********************************
    ## EOF
    ##---------------------------------
    #print(a3num)
    #print(a3num.shape)
    ##a3num = a3num - a3num.mean(axis=0)

    #a2dat = a3num.reshape(a3num.shape[0],-1)
    #print(a2dat.shape)
    #pca = PCA(n_components=4)
    #pca.fit(a2dat)


    #a2score = pca.transform(a2dat).T
    #print(a2score)
    #print(a2score.shape)

    #print(np.cumsum(pca.explained_variance_ratio_))

    #fig = plt.figure()
    #ax  = fig.add_subplot()
    #for i in range(a2score.shape[0]):
    #    plt.plot(lYear, a2score[i], label=i)

    #ax.legend()
    #plt.show()
    ##for itmp,dattmp in enumerate(pca.components_):
    ##    dattmp = dattmp.reshape(90,90)
    ##    draw_map(dattmp, stitle=itmp, cmbnd=np.arange(-0.02,0.02+0.001,0.005))

#%%
        ##***********************
        ## Figure: ensemble mean
        ##***********************
        #a3dat0 = ma.masked_less(a3dat0, 0)
        #a3dat1 = ma.masked_less(a3dat1, 0)

        #a2dif = (a3dat0.mean(axis=0) - a3dat1.mean(axis=0)) /len(lYear) *rp   # times/rp-year

        #a2tv, a2pv = scipy.stats.ttest_ind(a3dat1.filled(np.nan), a3dat0.filled(np.nan), axis=0, equal_var=False, nan_policy="omit")

        #a2fig = a2dif
        #a2fig = ma.masked_where(a2tcprec < 10, a2fig)

        #a2hatch = ma.masked_where(a2pv>0.05, a2fig)
        #a2hatch = ma.masked_where(a2tcprec<10, a2hatch)

        ##-- title, figure name -----
        #stitle = 'diff. tc-precip-extreme counts (times/%s-yr)\n'%(rp) + '%04d-%04d ens:%03d-%03d rad:%04dkm'%(lYear0[0],lYear1[-1],lens[0],lens[-1], radkm) 
        #figpath = figdir + '/map.timedif.tc-prec-extreme-counts.%04dkm.%s.rp-%03d.%04d-%04d.%s.png'%(radkm, region,rp,lYear0[0],lYear1[-1], scen)
        ##---------------------------
        #figmap   = plt.figure(figsize=(6,4))
        #axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
        #axmap.set_facecolor("0.8")
        ##-- grid lines ---
        #gl       = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        #xticks   = np.arange(-180, 180+1, 15)
        #yticks   = np.arange(-90,90, 15)
        #gl.xlocator = mticker.FixedLocator(xticks)
        #gl.ylocator = mticker.FixedLocator(yticks)
        #axmap.set_xticks(xticks, crs=ccrs.PlateCarree())
        #axmap.set_yticks(yticks, crs=ccrs.PlateCarree())
        ##-- set extent and coastlines----
        ##axmap.set_extent([lllon,urlon,lllat,urlat])
        #axmap.set_extent([105,150,10,45])
        #axmap.coastlines(color="k")

        ##-- color boundaries norm --------
        ##cmbnd = list(np.arange(-2.5,2.5+0.01,0.25))
        ##cmlabels = list(np.arange(-2.5,2.5+0.01,1))
        #if rp in [1,5]:
        #    cmbnd = list(np.arange(-0.25, 0.25+0.001, 0.05))
        #    cmlabels = [-0.20, -0.1, 0, 0.1, 0.2]
        #elif rp in [5]:
        #    cmbnd = list(np.arange(-0.3, 0.3+0.001, 0.05))
        #    cmlabels = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3]
        #elif rp in [10,20]:
        #    cmbnd = list(np.arange(-0.5, 0.5+0.001, 0.1))
        #    cmlabels = [-0.5,-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5]

        #mycm = 'RdBu_r'

        #cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        #cmaplist = [cmap(i) for i in range(cmap.N)]
        #cmap = matplotlib.colors.ListedColormap(cmaplist)
        #norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)

        ###-- contourf --------------
        ##X,Y = np.meshgrid(a1loncnt[x0:x1+1],a1latcnt[y0:y1+1])
        ###vmin, vmax = cmbnd[0], cmbnd[-1]
        ##im  = plt.contourf(X,Y,a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="both")

        #X,Y = np.meshgrid(a1lonbnd[x0:x1+1],a1latbnd[y0:y1+1])
        #vmin, vmax = cmbnd[0], cmbnd[-1]
        #im  = plt.pcolormesh(X,Y,a2fig, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
        #print(cmbnd)
        ##-- hatch --------------
        ##plt.contourf(X, Y, a2hatch, hatches=["////"], alpha=0.)
        #plt.pcolor(X, Y, a2hatch, hatch="///", alpha=0.)

        ##-- draw colorbar ------
        #cax = figmap.add_axes([0.84,0.2,0.02,0.6])
        #cbar= plt.colorbar(im, orientation='vertical', cax=cax)
        ##cbar.set_ticks(cbar.ax.get_yticks())
        #cbar.set_ticks(cmlabels)
        #cbar.set_ticklabels(cmlabels)

        ##-- Tiele -----
        #axmap.set_title(stitle)
        ##-- Save ------
        #plt.savefig(figpath)
        #plt.show()
        #print(figpath)

   
# %%
