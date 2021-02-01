# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
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
calcflag = True
#calcflag = False
figflag = True
figflag = False
iY = 1990
eY = 2010
#eY = 1995
lYear = range(iY,eY+1)
lMon  = range(1,12+1)
#lMon  = range(7,9+1)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB_NAT']
lscen   = ['HPB']
#lens    = list(range(1,50+1))
lens    = list(range(3,50+1))
noleap  = False

detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))

hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'


compbasedir= '/home/utsumi/temp/bams2020/composite'
figdir  = '/home/utsumi/temp/bams2020/tc-prec'
util.mk_dir(figdir)
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
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


const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = d4PDF.Lat()
const['Lon'] = d4PDF.Lon()

const['thsst']   = thsst + 273.15   # K
const['exrvort'] = exrvort
const['tcrvort'] = tcrvort 
const['thwcore'] = thwcore
const['thdura']  = thdura
const['thwind']  = thwind
const['thwdif']  = thwdif

a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")

#-- Load extreme criteria (NAT) ---
rp = 20  # return period, years
thdir = "/home/utsumi/temp/bams2020/extreme-prec-d/HPB_NAT.1990-2010"
thpath= thdir + "/prec.rp-020.WNP.npy"
a2th = np.load(thpath)    # original unit

#----------------------------------
for scen in lscen:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    if calcflag !=True: continue

    for ens in lens:
        #if (scen=='HPB')&(ens<4): continue # test

        print(('ens=',ens))
        wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    
        cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)


        for Year in lYear:
            #----------------------
            a2num = np.zeros([nyreg,nxreg], "int32")
            for Mon in lMon:
                print(scen,ens,Year,Mon)
                #-- test --------
                #if ((scen=="HPB")&(ens==1)&(datetime(Year,Mon,1,0)<=datetime(2000,9,1,0))):continue
                #if ((scen=="HPB_NAT")&(ens==2)&(datetime(Year,Mon,1,0)<=datetime(1991,8,1,0))):continue
                #----------------

                lDTime = util.ret_lDTime_fromYM([Year,Mon],[Year,Mon], timedelta(hours=24), hour0=0)

                #--- load precipitation ---
                a3prec_day = d4sfc.load_6hr_mon("PRECIPI", scen, ens, Year, Mon).reshape(-1,4,ny,nx).mean(axis=1)   # original unit
                a3prec_day = a3prec_day[:,y0:y1+1,x0:x1+1]
                #--- load TC dictionary ---
                YM = [Year,Mon]
                dexcxy, dtcxy  = cy.mkInstDictC_objTC(YM,YM,varname='vortlw')

                #--------------------------    

                for i,DTime in enumerate(lDTime):
                    #--------------------------
                    ltcxy = []
                    for dh in [0,6,12,18]:
                        ltcxy = ltcxy + dtcxy[DTime + timedelta(hours=dh)]

                    if len(ltcxy)==0: continue    

                    ltcx, ltcy, ltcvar = list(map(np.array,list(zip(*ltcxy))))

                    #-- Make mask ---
                    a2mask = detect_fsub.mk_a2mask_with_table(a2table.T, ltcx, ltcy, nx, ny).T
                    a2mask = a2mask[y0:y1+1,x0:x1+1]
                    a2mask = ma.masked_less_equal(a2mask,0)
                    a2mask = ma.masked_where(a3prec_day[i] <a2th, a2mask)
                    a2num  = a2num + a2mask.filled(0).astype("int32")

            #--- Save file ---------
            precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme"
            precdir = precbasedir + "/%s.%03d"%(scen,ens)
            util.mk_dir(precdir)

            numpath = precdir + "/num.%s.rp-%03d.%04d.npy"%(region, rp, Year)
            np.save(numpath, a2num.astype("int32"))
            print(numpath)

###*** Draw map *************
for scen in lscen:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    a3num = np.zeros([len(lens),nyreg,nxreg], "int32")
    for iens,ens in enumerate(lens):
        print(scen,ens)
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-extreme"

        a2numTmp = np.array([np.load( precbasedir + "/%s.%03d/num.%s.rp-%03d.%04d.npy"%(scen, ens, region, rp, Year)) for Year in lYear]).sum(axis=0)

        a3num[iens] = a2numTmp

    a2num = a3num.mean(axis=0) /len(lYear)*rp  # times/rp-year (e.g., rp=10 --> times/10-year)

    #-- Figure --------
    a2fig = a2num

    a1latreg = d4PDF.Lat()[y0:y1+1]
    a1lonreg = d4PDF.Lon()[x0:x1+1]
    X,Y = np.meshgrid(a1lonreg, a1latreg)
    BBox = [[0,100],[50,150]]
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

    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #--contourf ----------
    im = axmap.contourf(X,Y, a2fig, levels=cmbnd, cmap=mycm, extend="max")
    cbar = plt.colorbar(im)
    cbar.set_ticks(cmlabels)
    cbar.set_ticklabels(cmlabels)

    #-- title, figure name -----
    figpath = figdir + '/map.prec-tc.count-extreme.%s.%s.rp-%03d.%04d-%04d.png'%(region,scen,rp,iY,eY)
    stitle = 'Count TC extreme prec. times/%s-yr %s\n'%(rp, scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])

    axmap.set_title(stitle)


    plt.show()

    plt.savefig(figpath)
    print(figpath)


# %%
