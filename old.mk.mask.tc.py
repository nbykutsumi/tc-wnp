# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
#--------------------------------------
#calcbst= True
calcbst= False
#figbst = True
figbst = False

calcobj= True
#calcobj= False
figobj = True
#figobj = False

figmean = True

lens    = range(1,9+1)
#lens    = [1]

ctype = 'TC'
#ctype = 'ALL'


iDTime = datetime(2000,1,1,0)
#eDTime = datetime(2000,1,31,18)
#eDTime = datetime(2000,3,31,18)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2010,9,1,0)
#eDTime = datetime(2010,9,31,18)

dDTime = timedelta(hours=6)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]

cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = None
radkm = 500
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT']
lscen   = ['HPB']
#lscen   = ['HPB_NAT']
#lens    = range(1,9+1)
#lens    = [1]
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))


wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
figdir  = '/home/utsumi/temp/bams2020'
picklebasedir= '/home/utsumi/temp/bams2020/pickle'
d4pdfdir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
d4atm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir=d4pdfdir)

#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss =-9999
miss_int= -9999
dgrid = 5
a1latbnd  = np.arange(-90,90+0.01, dgrid)
a1lonbnd  = np.arange(0,360+0.01, dgrid)
ny = len(a1latbnd) - 1
nx = len(a1lonbnd) - 1
lonbnd0 = a1lonbnd[0]
latbnd0 = a1latbnd[0]

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]
[[lllat,lllon],[urlat,urlon]] = [[20,120],[47,150]]
a1latbndfig = np.arange(lllat, urlat+0.01, dgrid)
a1lonbndfig = np.arange(lllon, urlon+0.01, dgrid)

ixregion= bisect_left(a1lonin, lllon)
exregion= bisect_left(a1lonin, urlon)
iyregion= bisect_left(a1latin, lllat)
eyregion= bisect_left(a1latin, urlat)


#nyfig = a1latbnd.shape[0] -1
#nxfig = a1lonbnd.shape[0] -1

#**************************

##**************************
##Read best track
##---------------------------
##iDTimeBst = datetime(1980,1,1,0)
##eDTimeBst = datetime(2019,12,31)
#iDTimeBst = datetime(2000,1,1,0)
#eDTimeBst = datetime(2010,12,31)
#
#
#lDTimeBst = util.ret_lDTime(iDTimeBst,eDTimeBst,timedelta(hours=6))
#if calcbst is True:
#    dbst   = {}
#    bst    = IBTrACS.IBTrACS()
#    dlonlat = bst.ret_dlonlat(iDTimeBst,eDTimeBst)
#    llonlat = []
#    for lonlat in dlonlat.values():
#        llonlat = llonlat + lonlat
#    #-- Map --------
#    a2count = np.zeros([ny,nx],int32) 
#    for (lon,lat) in llonlat:
#        ix = int((lon-lonbnd0)/dgrid)
#        iy = int((lat-latbnd0)/dgrid)
#        if (ix<0)or(ix>nx-1)or(iy<0)or(iy>ny-1): continue
#        a2count[iy,ix] = a2count[iy,ix] + 1
#
#    a2freq = a2count.astype('float32')/len(lDTimeBst)
#
#    dbst['a2count'] = a2count
#    dbst['a2freq' ] = a2freq
#    dbst['iDTime']= iDTimeBst
#    dbst['eDTime']= eDTimeBst
#    dbst['a1latbnd'] = a1latbnd
#    dbst['a1lonbnd'] = a1lonbnd
#    dbst['dgrid']    = dgrid
#
##-- Save --------
#bstPath= outDir + '/freq.tc.bst.pickle'
#if calcbst is True:
#    with open(bstPath,'wb') as f:
#        pickle.dump(dbst, f)
#    print bstPath
#
##-- Figure (best track)
#if figbst is True:
#    with open(bstPath, 'rb') as f:
#        dbst = pickle.load(f)
#
#    a2fig = dbst['a2freq']
#    a2fig = ma.masked_less(a2fig,0)*4*365  # times/year
#    dpara = {}
#    dpara['title'] = 'Prob. of existence (Best track)'
#    figdir  = '/home/utsumi/temp/bams2020'
#    dpara['figpath'] = figdir + '/map.freq.tc.bst.png'
#    dpara['cmbnd']   = cmbnd
#
#    draw_map(a2fig, dpara)


#************************************
# d4PDF (Objective detection)
#************************************
#lexrvort= 3.7*1.0e-5 * np.array([1, 1.3])
#lthwcore= [0.2]

#lexrvort= np.array([2.0, 2.5, 3.0])*1.0e-5
#lthwcore= [0.2]

#lthsst  = [25,27,29]
lthsst  = [27]
#lthsst  = [15,25]
#lthsst  = [-273.15]
#lexrvort= np.array([3.7])*1.0e-5
#lexrvort= np.array([2.5,3.0])*1.0e-5
#lexrvort= np.array([2.0, 2.5])*1.0e-5
#lexrvort= np.array([-9999])*1.0e-5
#lexrvort= np.array([6.0])*1.0e-5
lexrvort = np.array([3])*1.0e-5
ltcrvort = np.array([4.8])*1.0e-5
#lthwcore= [0.0]
lthwcore= [0.3]
#lthwcore= [0, 0.5]
#lthwcore= [-99999]
#lthwcore= [-99]
#lthdura = [48,60,72]
lthdura = [48]
#lthdura = [6,12,24,36]
#lthdura = [-99999,-6,12,24,36]

lkey = [[thsst,exrvort,tcrvort,thwcore,thdura]
        for thsst   in lthsst
        for exrvort in lexrvort
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thdura in lthdura
        ]

for (thsst,exrvort,tcrvort,thwcore,thdura) in lkey:
    const  = ConstCyclone.Const(prj=prj, model=model)
    const['Lat'] = d4PDF.Lat()
    const['Lon'] = d4PDF.Lon()

    const['thsst']   = thsst + 273.15   # K
    const['exrvort'] = exrvort
    const['tcrvort'] = tcrvort 
    const['thwcore'] = thwcore
    const['thdura']  = thdura
    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    for scen in lscen:
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thdura)
        #slabel = 'wceach-%s-%s-%s-sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(ctype, expr, scen, thsst, exrvortout, tcrvortout, thwcore, thdura)

        for ens in lens:
            if calcobj is True:
                print 'ens=',ens
                lprec = []
                lwind = []
                #-------------------
                wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
    
                cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)

                ym0 = (-9999,-99)
                for DTime in lDTime:
                    y0,m0 = DTime.timetuple()[:2]

                    #--- load TC dictionary ---
                    if (y0,m0) != ym0:
                        ym0 = (y0,m0)
                        dexcxy, dtcxy  = cy.mkInstDictC_objTC(ym0,ym0,varname='vortlw')

                    #--------------------------
                    ltcxy = dtcxy[DTime]

                    if len(ltcxy)==0: continue
                    ltcx, ltcy, ltcvar = map(list,zip(*ltcxy))

                    #-- projection on map -----
                    a2loc = np.full([nyin,nxin], miss, 'float32')
                    a2loc[ltcy,ltcx] = 1

                    a2flag = np.zeros([nyin,nxin], 'float32')
                    a2flag[iyregion:eyregion+1,ixregion:exregion+1] = 1
                    a2loc = ma.masked_where(a2flag !=1, a2loc).filled(miss)

                    a2tcmask = detect_fsub.mk_territory(a2loc.T, a1lonin, a1latin, radkm*1000., imiss=miss, omiss=miss).T

                    if a2tcmask.max()<1: continue

                    a2prec = d4sfc.load_6hr(vname='PRECIPI', scen=scen, ens=ens, DTime=DTime) *60*60
                    a2u    = d4sfc.load_6hr(vname='UAOPN', scen=scen, ens=ens, DTime=DTime)
                    a2v    = d4sfc.load_6hr(vname='VAOPN', scen=scen, ens=ens, DTime=DTime)
                    a2wind = np.sqrt(a2u**2 + a2v**2)

                    #-- screen small variables --
                    a2maskrad  = ma.masked_less(a2tcmask,0).mask

                    a2prec = ma.masked_where(a2maskrad, a2prec)

                    a2wind = ma.masked_where(a2maskrad, a2wind)

                    lprec = lprec + a2prec.compressed().tolist()
                    lwind = lwind + a2wind.compressed().tolist()
                    #print len(lprec)
                    print scen, ens, DTime
                #--- Save file ---------
                pickledir = picklebasedir + '/rad-%d-%s'%(radkm, slabel)
                util.mk_dir(pickledir)

                picklepath =  pickledir + '/prec.%s.%03d.pickle'%(scen, ens)
                with open(picklepath,'wb') as f:
                    pickle.dump(lprec, f)

                picklepath =  pickledir + '/wind.%s.%03d.pickle'%(scen, ens)
                with open(picklepath,'wb') as f:
                    pickle.dump(lwind, f)

                print picklepath
    # %%
    #---- Draw PDF ---------------
    dprec = {}
    dwind = {}
    for scen in lscen:
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thdura)
        pickledir = picklebasedir + '/rad-%d-%s'%(radkm, slabel)

        #for ens in lens:
        for ens in [1]:
            if figobj is not True: continue
            picklepath =  pickledir + '/prec.%s.%03d.pickle'%(scen, ens)
            with open(picklepath, 'rb') as f:
                dprec[scen,ens] = pickle.load(f)

            picklepath =  pickledir + '/wind.%s.%03d.pickle'%(scen, ens)
            with open(picklepath, 'rb') as f:
                dwind[scen,ens] = pickle.load(f)
                


    #---- PDF for Prec ----
    for scen in lscen:
        for ens in lens:

            a1prec = ma.masked_less_equal(dprec[scen,ens]*60*60,1).compressed()
            print len(a1prec)
            plt.hist(a1prec)
            figpath = '/home/utsumi/temp/bams2020/temp.png'
            plt.savefig(figpath)
            print figpath
            sys.exit()
        
         
