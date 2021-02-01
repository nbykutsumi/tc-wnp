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
figflag  = True
#figflag  = False
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
#lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
lscen   = ['HPB_NAT']
lens    = list(range(1,50+1))
#lens    = list(range(1,1+1))
#wday   = 1  # 1, 3
wday   = 3  # 1, 3
ldh    = {1:[0,6,12,18], 3:np.arange(-24,42+1,6)}[wday]
print(len(ldh))
region = "WNP"
dBBox = {"WNP":[[0,100],[50,150]]}
[[lllat,lllon],[urlat,urlon]] = dBBox[region]

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
figdir  = '/home/utsumi/temp/bams2020/fig/map-tc-count-perday'
util.mk_dir(figdir)
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
#a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', mipyss_fill=-9999.)
#----------------------------------
a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

ny, nx = 320,640
miss =-9999
miss_int= -9999


y0 = bisect_left(a1lat, lllat)
y1 = bisect_left(a1lat, urlat)
x0 = bisect_left(a1lon, lllon)
x1 = bisect_left(a1lon, urlon)
nyreg = y1-y0+1
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

#radkm = 200
radkm = 500

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

if radkm == 500:
    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")
elif radkm == 200:
    a2table = np.load("./tab.dydx4mask.d4PDF.nrady-003.0200km.npy")
else:
    print("check radkm",radkm)
    sys.exit()

for scen in lscen:
    if calcflag !=True: continue

    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    for ens in lens:
        #if (scen=='HPB')&(ens<4): continue # test

        print(('ens=',ens))
        wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    
        cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)


        for Year in lYear:
            ndays = (datetime(Year,12,31,0) - datetime(Year,1,1,0)).days +1

            iday = -1 
            a3num = np.full([ndays,nyreg,nxreg], miss, "int8")

            for Mon in lMon:
                print(scen,ens,Year,Mon)
                #-- test --------
                #if ((scen=="HPB")&(ens==1)&(datetime(Year,Mon,1,0)<=datetime(2000,9,1,0))):continue
                #if ((scen=="HPB_NAT")&(ens==2)&(datetime(Year,Mon,1,0)<=datetime(1991,8,1,0))):continue
                #----------------

                lDTime = util.ret_lDTime_fromYM([Year,Mon],[Year,Mon], timedelta(hours=24), hour0=0)

                #--- load TC dictionary ---
                if wday ==1:
                    iYM = [Year,Mon]
                    eYM = [Year,Mon]
                elif wday==3:
                    iYM = util.shift_YM(Year,Mon,-1)
                    eYM = util.shift_YM(Year,Mon,+1)

                try:
                    dexcxy, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
                except:
                    dexcxy, dtcxy  = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname='vortlw')

                #--------------------------    
                for i,DTime in enumerate(lDTime):
                    iday = iday+1
                    #--------------------------

                    a2num = np.zeros([nyreg,nxreg],"int8")
                    for dh in ldh:
                        try:
                            ltcxy = dtcxy[DTime + timedelta(hours=int(dh))]
                        except KeyError:
                            continue

                        if len(ltcxy)==0: continue    
    
                        ltcx, ltcy, ltcvar = list(map(np.array,list(zip(*ltcxy))))
    
                        #-- Make mask ---
                        a2mask = detect_fsub.mk_a2mask_with_table(a2table.T, ltcx, ltcy, nx, ny).T
                        a2num = a2num + a2mask[y0:y1+1,x0:x1+1].astype("int8")

                    a3num[iday] = a2num
            #--- Save file ---------
            precbasedir = "/home/utsumi/temp/bams2020/tc-prec-%04dkm"%(radkm)
            precdir = precbasedir + "/%s.%03d"%(scen,ens)
            util.mk_dir(precdir)

            countpath = precdir + "/timesper-%dday.%04d.npy"%(wday,Year)

            np.save(countpath, a3num.astype("int8"))

            print(countpath)


##*** Draw map *************
for scen in lscen:
    if figflag !=True: continue

    a3ave = np.full([len(lens),nyreg,nxreg], -9999).astype("float64")
    for iens,ens in enumerate(lens):
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-%04dkm"%(radkm)
        precdir = precbasedir + "/%s.%03d"%(scen,ens)

        a3pass = np.full([len(lYear),nyreg,nxreg], -9999)
        a3even = np.full([len(lYear),nyreg,nxreg], -9999)

        for i,Year in enumerate(lYear):
            print(ens,Year)
            a3count = ma.masked_equal(np.load(precdir + "/timesper-%dday.%04d.npy"%(wday,Year)), 0)

            a2passTmp = a3count.sum(axis=0)
            a2evenTmp = a3count.count(axis=0)

            a3pass[i] = a2passTmp
            a3even[i] = a2evenTmp

        a2pass   = a3pass.sum(axis=0)
        a2even   = a3even.sum(axis=0)
        a2aveTmp = (ma.masked_where(a2even==0, a2pass.astype("float64")) / a2even).filled(np.nan)


        a3ave[iens] = a2aveTmp

    a2fig = ma.masked_invalid(a3ave).mean(axis=0)



    [[lllat,lllon],[urlat,urlon]] = dBBox[region]
    a1lat = d4PDF.Lat()
    a1lon = d4PDF.Lon()
    y0 = bisect_left(a1lat, lllat)
    y1 = bisect_left(a1lat, urlat)
    x0 = bisect_left(a1lon, lllon)
    x1 = bisect_left(a1lon, urlon)
   
    X,Y = np.meshgrid(a1lon[x0:x1+1], a1lat[y0:y1+1])

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
    cmbnd = list(np.arange(0,3+0.01,0.1))
    #cmlabels = list(np.arange(0,50+1,10))

    cmap   = plt.cm.get_cmap(mycm, len(cmbnd)+1)  # define the colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = matplotlib.colors.ListedColormap(cmaplist)
    norm = matplotlib.colors.BoundaryNorm(cmbnd, ncolors=cmap.N, clip=False)
    #--extent and coast lines ------
    axmap.set_extent([lllon,urlon,lllat,urlat])
    axmap.coastlines()

    #--contourf ----------
    im = axmap.contourf(X,Y, a2fig, levels=cmbnd, cmap=mycm, extend="max")
    plt.colorbar(im)

    #-- title, figure name -----
    figpath = figdir + '/map.count-per%dday.%s.%04d-%04d.png'%(wday,scen,iY,eY)
    stitle = 'average count/%dday %s\n'%(wday,scen) + '%04d-%04d ens:%03d-%03d'%(iY,eY,lens[0],lens[-1])

    axmap.set_title(stitle)


    plt.show()

    plt.savefig(figpath)
    print(figpath)

# %%
