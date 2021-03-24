# %%
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
#%matplotlib inline

import cartopy.crs as ccrs
import matplotlib.ticker as mticker

import numpy as np
import scipy
import d4PDF
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from bisect import bisect_left
from importlib import import_module
from numpy import ma

calcflag = True
#calcflag = False
figflag = True

prj     = "d4PDF"
model   = "__"
expr    = 'XX'
lscen   = ['HPB']  # run={expr}-{scen}-{ens}
#lscen   = ['HPB_NAT'] # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT'] # run={expr}-{scen}-{ens}

#lens    = list(range(1,50+1))
#lens    = list(range(1,20+1))
lens    = [1]
miss_out = -9999.

#iY,eY = 1990,2010
iY,eY = 2001,2010
#iY,eY = 1990,1990
lY  = range(iY,eY+1)
lM  = range(1,12+1)
ny,nx = 320, 640
hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'

detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)

a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

region= "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

y0 = bisect_left(a1lat, lllat)
y1 = bisect_left(a1lat, urlat)
x0 = bisect_left(a1lon, lllon)
x1 = bisect_left(a1lon, urlon)

nyreg = y1-y0+1    # 
nxreg = x1-x0+1

print(nxreg,nyreg)
wy = 30
ly00 = np.arange(y0,y1+1)[::wy]

ndays = (datetime(eY,12,31)-datetime(iY,1,1)).days +1
lrp = [1,5,10,20]
dpercent = {rp:100-100/(365*rp) for rp in lrp}
print(dpercent)
for scen in lscen:
    if calcflag != True: continue

    for y00 in ly00:
        a3dat = np.zeros([ndays*len(lens), wy, nxreg], "float32")
        print(a3dat.shape)
        for iens, ens in enumerate(lens):
            for Year in lY:
                for Mon in lM:

                    a3temp = d4sfc.load_6hr_mon("PRECIPI", scen, ens, Year, Mon).reshape(-1,4,ny,nx).mean(axis=1)[:,y00:y00+wy,x0:x1+1] 
                    print(scen,ens,Year,Mon,a3temp.shape, ndays)

                    iday = (datetime(Year,Mon,1) - datetime(lY[0],1,1)).days
                    a3dat[iens*ndays+iday:iens*ndays+iday+a3temp.shape[0],:,:] = a3temp

        for rp in lrp:
            percent = dpercent[rp]
            a2ptile = np.percentile(a3dat, percent, axis=0)
        
            #--- Save -----
            exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/%s.ens-%03d-%03d-yr-%04d-%04d"%(scen,lens[0],lens[-1],iY,eY)
            util.mk_dir(exdir)
            expath= exdir + "/prec.rp-%03d.%s.y-%03d-%03d.npy"%(rp, region, y00, y00+wy-1)
            np.save(expath, a2ptile)
            print(expath)
#-- Join ------------
for scen in lscen:
    for rp in lrp:
        #exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/%s.%04d-%04d"%(scen,iY,eY)
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/%s.ens-%03d-%03d-yr-%04d-%04d"%(scen,lens[0],lens[-1],iY,eY)

        a2exp = np.concatenate([np.load(exdir + "/prec.rp-%03d.%s.y-%03d-%03d.npy"%(rp, region, y00, y00+wy-1)) for y00 in ly00], axis=0)

        joinpath = exdir + "/prec.rp-%03d.%s.npy"%(rp, region)
        np.save(joinpath, a2exp)
        print(joinpath)

#***** Figure *********
for scen in lscen:
    if figflag != True: continue
    for rp in lrp:
        exdir = "/home/utsumi/temp/bams2020/extreme-prec-d/%s.ens-%03d-%03d-yr-%04d-%04d"%(scen,lens[0],lens[-1],iY,eY)
        joinpath = exdir + "/prec.rp-%03d.%s.npy"%(rp, region)
        a2fig = np.load(joinpath) * 60*60*24  # mm/sec --> mm/day

        fig = plt.figure(figsize=(6,4))
        axmap  = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())

        gl = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
        xticks = np.arange(-180,180+1,15)
        yticks = np.arange(-90,90+1,15)
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)
        axmap.set_xticks(xticks, crs = ccrs.PlateCarree())
        axmap.set_yticks(yticks, crs = ccrs.PlateCarree())

        #-- Make new colormap (white at the lower end) --
        upper = matplotlib.cm.jet(np.arange(256))
        lower = np.ones((int(256/4),4))
        for i in range(3):
          lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
        mycm = np.vstack(( lower, upper ))
        mycm = matplotlib.colors.ListedColormap(mycm, name='myColorMap', N=mycm.shape[0])

        #-- color boundaries norm ------
        #clevs = range(-20,20+1,4)
        cmbnd = [0,40,80,120,160,200,300,400,500]


        cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = matplotlib.colors.ListedColormap(cmaplist)
        norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
        #--extent and coast lines ------
        dbboxfig = {"WNP":[[0,100],[50,150]]}
        [[lllatfig,lllonfig],[urlatfig,urlonfig]] = dbboxfig[region]
        axmap.set_extent([lllonfig,urlonfig,lllatfig,urlatfig])
        axmap.coastlines()

        #-- title, figure name -----
        figdir  = "/home/utsumi/temp/bams2020/fig/map-prec"
        figpath = figdir + "/map.prec.%s.rp-%03d.ens-%03d-%03d-yr-%04d-%04d.%s.png"%(scen,rp, lens[0],lens[-1],iY,eY,region)
        stitle = '%s-year R.P. mm/day %s\n'%(rp, scen) + 'ens:%03d-%03d yr:%04d-%04d'%(lens[0],lens[-1],iY,eY)

        axmap.set_title(stitle)

        #-- contourf -------------------
        a1latreg = a1lat[y0:y1+1]
        a1lonreg = a1lon[x0:x1+1]

        X,Y = np.meshgrid(a1lonreg, a1latreg)
        im = axmap.contourf(X,Y, a2fig, cmap=cmap, norm=norm, levels=cmbnd, extend="max")
        plt.colorbar(im)
        plt.show()

        plt.savefig(figpath)
        print(figpath)



    




# %%

# %%
