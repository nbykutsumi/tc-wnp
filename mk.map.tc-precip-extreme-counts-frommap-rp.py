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
#--------------------------------------
#lrp = [1,5,10,20]  # return period, years
lrp = [1]  # return period, years
radkm  = 500 # 200, 500

calcflag = True
#calcflag = False
#figflag = True
figflag = False
#iY = 1990
#eY = 2010

iY = 2001
eY = 2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT']
lscen   = ['HPB_NAT']
#lscen   = ['HPB']
#lens    = list(range(1,50+1))
#lens    = list(range(1,20+1))
lens    = list(range(1,1+1))

#iYbase = 1990
#eYbase = 2010
#lensbase = list(range(1,50+1))
#scenbase = "HPB_NAT"

iYbase = 2001
eYbase = 2010
#lensbase = list(range(1,20+1))
lensbase = list(range(1,1+1))
#scenbase = "HPB"
scenbase = "HPB_NAT"

detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
    srcbasedir = "/tank/utsumi/hometemp/bams2020"
    
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'
    srcbasedir = "/home/utsumi/mnt/lab_tank/utsumi/hometemp/bams2020"

compbasedir= '/home/utsumi/temp/bams2020/composite'
figdir  = '/home/utsumi/temp/bams2020/tc-prec'
util.mk_dir(figdir)
#d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640

region= "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

y0 = bisect_left(a1lat, lllat)
y1 = bisect_left(a1lat, urlat)
x0 = bisect_left(a1lon, lllon)
x1 = bisect_left(a1lon, urlon)

nyreg = y1-y0+1    #
nxreg = x1-x0+1

#************************************
# d4PDF (Objective detection)
#************************************
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

if radkm == 500:
    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")
elif radkm == 200:
    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-003.0200km.npy")

#----------------------------------
for rp in lrp:
    #-- Load extreme criteria (NAT) ---
    #thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/HPB_NAT.1990-2010"
    #thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/HPB_NAT.ens-001-050-yr-1990-2010"
    thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/%s.ens-%03d-%03d-yr-%04d-%04d"%(scenbase, lensbase[0], lensbase[-1], iYbase, eYbase)

    thpath= thdir + "/prec.rp-%03d.WNP.npy"%(rp)
    a2th = np.load(thpath)    # original unit



    for scen in lscen:
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    
        if calcflag !=True: continue
    
        for ens in lens:
            #if (scen=='HPB')&(ens<4): continue # test
    
            print(('ens=',ens))
            #wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
        
            #cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)
    
    
            for Year in lYear:
    
                srcdir = srcbasedir + "/tc-prec-%04dkm/%s.%03d"%(radkm,scen,ens)
                precpath= srcdir + "/prec-tc.%s.%04d.npy"%(region,Year)
                a3tcprec_day = np.load(precpath)
                a2num = ma.masked_greater_equal(a3tcprec_day, np.broadcast_to(a2th, a3tcprec_day.shape)).mask.sum(axis=0).astype("int32")
    
                #--- Save file ---------
                #precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
                precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm-base-%s-ens-%03d-%03d-yr-%04d-%04d"%(radkm, scenbase, lensbase[0], lensbase[-1], iYbase, eYbase)

                precdir = precbasedir + "/%s.%03d"%(scen,ens)
                util.mk_dir(precdir)
    
                numpath = precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)
                np.save(numpath, a2num.astype("int32"))
                print(numpath)

    ###*** Draw map *************
    for scen in lscen:
        if figflag != True: continue
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    
        a3num = np.zeros([len(lens),nyreg,nxreg], "int32")
        for iens,ens in enumerate(lens):
            print(scen,ens)
            #precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm"%(radkm)
            precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme-%04dkm-base-%s-ens-%03d-%03d-yr-%04d-%04d"%(radkm, scenbase, lensbase[0], lensbase[-1], iYbase, eYbase)
    
            a2numTmp = np.array([np.load( precbasedir + "/%s.%03d/num.%s.rp-%03d.%04d.npy"%(scen, ens, region, rp, Year)) for Year in lYear]).sum(axis=0)
    
            a3num[iens] = a2numTmp
    
        a2num = a3num.mean(axis=0) /len(lYear)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)
    
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
        if rp in [1,5,10,20]:
            cmbnd = list(np.arange(0,1+0.001, 0.1))
            cmlabels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0]

            #cmbnd = list(np.arange(0,0.9+0.001, 0.1))
            #cmlabels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


            #cmbnd = list(np.arange(0,1+0.001, 0.2))
            #cmlabels = [0, 0.2, 0.4, 0.6, 0.8, 1.0]



        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
        #--extent and coast lines ------
        axmap.set_extent([lllon,urlon,lllat,urlat])
        axmap.coastlines()
    
        ##--contourf ----------
        #a1latreg = d4PDF.Lat()[y0:y1+1]
        #a1lonreg = d4PDF.Lon()[x0:x1+1]
        #X,Y = np.meshgrid(a1lonreg, a1latreg)
        #im = axmap.contourf(X,Y, a2fig, levels=cmbnd, cmap=mycm, extend="max")
        #cbar = plt.colorbar(im)
        #cbar.set_ticks(cmlabels)
        #cbar.set_ticklabels(cmlabels)

        #-- meshgrid ---------
        a1latbndreg = d4PDF.LatBnd()[y0:y1+1]
        a1lonbndreg = d4PDF.LonBnd()[x0:x1+1]
        X,Y = np.meshgrid(a1lonbndreg, a1latbndreg)
        im = axmap.pcolormesh(X,Y, a2fig, cmap=mycm, norm=norm)
        cbar = plt.colorbar(im)
        cbar.set_ticks(cmlabels)
        cbar.set_ticklabels(cmlabels)

        #-- title, figure name -----
        figpath = figdir + '/map.prec-tc.count-extreme.%04dkm.base-%s-ens-%03d-%03d-yr-%04d-%04d.%s.%s.rp-%03d.%04d-%04d.png'%(radkm,scenbase, lensbase[0], lensbase[-1], iYbase, eYbase,region,scen,rp,iY,eY)

        stitle = 'Count TC extreme prec. times/%s-yr %s\n'%(rp, scen) + '%04d-%04d ens:%03d-%03d %04dkm'%(iY,eY,lens[0],lens[-1], radkm) + "\nbase-%s b-ens:%03d-%03d b-yr-%04d-%04d"%(scenbase, lensbase[0], lensbase[-1], iYbase, eYbase)
    
        axmap.set_title(stitle)
    
    
        plt.show()
    
        plt.savefig(figpath)
        print(figpath)


# %%
