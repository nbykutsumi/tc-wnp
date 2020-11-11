# %%
import matplotlib
#%matplotlib inline
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import calendar
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
import IBTrACS
import Cyclone
import ConstCyclone
import d4PDF
#--------------------------------------
calcobs = False
#calcobs = True

calcsim = False
#calcsim = True

figmon = True
#figmon = False
figyear= True
#figyear= False

#target='track'
target='point'

ly_sim = range(2000,2010+1)
#ly_obs = range(1989,1989+1)
#ly_obs = range(1980,2018+1)
ly_obs = range(1980,2018+1)
lmon= range(1,12+1)

[[lllat,lllon],[urlat,urlon]] = [[20,120],[47,150]]
########[[lllat,lllon],[urlat,urlon]] = [[30,120],[47,150]]
#[[lllat,lllon],[urlat,urlon]] = [[0,120],[47,150]]

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
lens    = range(1,9+1)
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)
#countbasedir = '/home/utsumi/temp/bams2020/count-%s/%s'%(target,sdomain)
countbasedir = '/home/utsumi/temp/bams2020/count-point/%s'%(sdomain)
figdir  = '/home/utsumi/temp/bams2020/fig/%s'%(sdomain)
util.mk_dir(figdir)
#**********************************
# Obs
#----------------------------------
ib = IBTrACS.IBTrACS()
for Year in ly_obs:
    if calcobs is False: continue
    for Mon in lmon:
        print Year,Mon
        eday = calendar.monthrange(Year,Mon)[1]
        idtime=datetime(Year,Mon,1,0)
        edtime=datetime(Year,Mon,eday,23)
        lbst = ib.ret_dlonlat(idtime,edtime,lvaridx=[0]).values()
        lbst = [item for l in lbst for item in l]  # flatten
        a1lon, a1lat, a1tcid = map(np.array, zip(*lbst))
        #-- Screening region ----
        a1flaglat = ma.masked_inside(a1lat, lllat, urlat).mask
        a1flaglon = ma.masked_inside(a1lon, lllon, urlon).mask
        a1flag = a1flaglat * a1flaglon

        if a1flag is np.bool_(False):
            stcid = np.array([])
        else:
            #------------------------
            a1tcid= a1tcid[a1flag]
            #if target=='track':
            #    stcid = set(a1tcid)
            #    stcid = np.sort(np.array(list(stcid)))
            #elif target=='point':
            #    stcid = np.sort(a1tcid)
            #else:
            #    print 'check target',target
            #    sys.exit()
            stcid = np.sort(a1tcid)

        countdir = countbasedir + '/obs'
        util.mk_dir(countdir)
        obspath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)
        np.save(obspath, stcid)
        print 'OBS',Year,Mon,len(stcid)
        print obspath

#sys.exit()
#----------------------------------
miss =-9999
miss_int= -9999

a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

lthsst = [27, 28]
exrvort= 3.0*1e-5
#ltcrvort= np.array([5, 6, 7])*1e-5
ltcrvort= np.array([5])*1e-5
lthwcore= [3, 4, 5]
#lthwcore= [3]
thdura = 36
lthwind = [15]
lthwdif = [-9999., -30]
#lthwdif = [-20]
#lthwdif = [-10.]
#lthwdif = [-15.]

lkey = [[thsst, tcrvort,thwcore,thwind,thwdif]
        for thsst   in lthsst
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thwind  in lthwind
        for thwdif  in lthwdif
        ]

for (thsst,tcrvort,thwcore,thwind,thwdif) in lkey:

    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    const  = ConstCyclone.Const(prj=prj, model=model)
    const['Lat']     = a1lat
    const['Lon']     = a1lon
    const['thsst']   = thsst + 273.15   # K
    const['exrvort'] = exrvort
    const['tcrvort'] = tcrvort
    const['thwcore'] = thwcore
    const['thwind']  = thwind
    const['thwdif']  = thwdif
    const['thdura']  = thdura


    d2_sim = {}
    d2hist_sim = {}
    for scen in lscen:
        d2hist_sim[scen] = []
        for ens in lens:
            if calcsim is False: continue

            wsDir = wsbaseDir + '/%s-%s-%03d'%(expr,scen,ens)
            cy = Cyclone.Cyclone(baseDir=wsDir, const=const)

            for Year in ly_sim:
                for Mon in lmon:
                    #--- Load file ---------
                    _, dxyipos   = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname='ipos')
                    _, dxyidate  = cy.mkInstDictC_objTC([Year,Mon],[Year,Mon],varname='idate')

                    lxyipos = [item for l in dxyipos.values() for item in l]  # flatten
                    lxyidate = [item for l in dxyidate.values() for item in l]  # flatten

                    lx,ly, lipos = zip(*lxyipos)
                    _,  _, lidate= zip(*lxyidate)

                    lx = np.array(lx)
                    ly = np.array(ly)

                    a1latTmp = a1lat[ly]
                    a1lonTmp = a1lon[lx]

                    #-- Region mask ---------
                    a1flaglat = ma.masked_inside(a1latTmp, lllat, urlat).mask
                    a1flaglon = ma.masked_inside(a1lonTmp, lllon, urlon).mask

                    a1flag = a1flaglat* a1flaglon

                    if a1flag is np.bool_(False):
                        stcid_sim = np.array([])

                    else:
                        a1idate = np.array(lidate)[a1flag]
                        a1ipos  = np.array(lipos )[a1flag]

                        a1tcid_sim = ['%s-%s'%(x,y) for (x,y) in zip(a1idate, a1ipos)]

                        #if target=='track':
                        #    stcid_sim = set(a1tcid_sim)
                        #    stcid_sim = np.sort(np.array(list(stcid_sim)))
                        #elif target=='point':
                        #    stcid_sim = np.sort(a1tcid_sim)
                        #else:
                        #    print 'check target',target
                        #    sys.exit()

                        stcid_sim = np.sort(a1tcid_sim)

                    countdir = countbasedir + '/%s/%s.%03d'%(slabel, scen, ens)
                    simpath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)
                    util.mk_dir(countdir)

                    np.save(simpath, stcid_sim)
                    print simpath

    #********************************
    # Monthly count
    #--------------------------------
    # Obs (monthly count)
    #-------------
    a2num_obs = []
    for Year in ly_obs:
        a1num_tmp = []
        for Mon in lmon:
            countdir = countbasedir + '/obs'
            util.mk_dir(countdir)
            obspath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)

            a1tcid_tmp = np.load(obspath)
            if target=='track':
                a1tcid_tmp = set(a1tcid_tmp)

            a1num_tmp.append(len(a1tcid_tmp))

        a2num_obs.append(a1num_tmp)

    if target=='point':
        nsteps = 30*4
        a2num_obs = np.array(a2num_obs) / float(nsteps)

    a2num_obs = np.array(a2num_obs)
    #-------------
    # Sim (monthly count)
    #-------------
    d3num_sim = {}
    for scen in lscen:
        a3num_sim = []
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

        for ens in lens:
            for Year in ly_sim:
                for Mon in lmon:
                    countdir = countbasedir + '/%s/%s.%03d'%(slabel, scen, ens)
                    simpath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)

                    try:
                        a1tcid_tmp = np.load(simpath)
                    except IOError:
                        a1tcid_tmp = np.array([])

                    if target=='track':
                        a1tcid_tmp = set(a1tcid_tmp)

                    a3num_sim.append(len(a1tcid_tmp))
        a3num_sim = np.array(a3num_sim).reshape(len(lens),len(ly_sim),-1)

        if target=='point':
            nsteps = 30*4
            a3num_sim = a3num_sim / float(nsteps)

        d3num_sim[scen] = a3num_sim

    #---------------------------
    # Figure (mounthly count)
    #---------------------------
    if figmon == True:
        cm_ens = {'HPB':'0.5','HPB_NAT':'blue'}
        cm_ave = {'HPB':'tomato',  'HPB_NAT':'steelblue'}

        fig = plt.figure(figsize=(7,4))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        a1x = np.arange(1,12+1)

        a1num_obs_mean= a2num_obs.mean(axis=0)
        a1num_his_mean= d3num_sim['HPB'].mean(axis=(0,1))
        a1num_nat_mean= d3num_sim['HPB_NAT'].mean(axis=(0,1))

        yerr_obs_lw =-np.percentile(a2num_obs, 25, axis=0) + a1num_obs_mean
        yerr_obs_up = np.percentile(a2num_obs, 75, axis=0) - a1num_obs_mean

        yerr_his_lw =-np.percentile(d3num_sim['HPB'], 25, axis=(0,1)) + a1num_his_mean

        yerr_his_up = np.percentile(d3num_sim['HPB'], 75, axis=(0,1)) - a1num_his_mean
        yerr_nat_lw =-np.percentile(d3num_sim['HPB_NAT'], 25, axis=(0,1)) + a1num_nat_mean
        yerr_nat_up = np.percentile(d3num_sim['HPB_NAT'], 75, axis=(0,1)) - a1num_nat_mean


        wbar = 0.25*(a1x[1] - a1x[0])
        a1x_obs = a1x - wbar
        a1x_nat = a1x 
        a1x_his = a1x + wbar

        ax.bar(a1x_obs, a1num_obs_mean, width=wbar, yerr=[yerr_obs_lw, yerr_obs_up], color='0.3', ecolor='k')
        ax.bar(a1x_nat, a1num_nat_mean, width=wbar, yerr=[yerr_nat_lw, yerr_nat_up], color='steelblue', ecolor='b')
        ax.bar(a1x_his, a1num_his_mean, width=wbar, yerr=[yerr_his_lw, yerr_his_up], color='tomato', ecolor='r')

        stitle = target + '\n' + slabel + '\n' +sdomain
        plt.title(stitle)
        #plt.show()
        figpath = figdir + '/count.mon.%s.%s.png'%(target, slabel)
        plt.savefig(figpath)
        print figpath


    #********************************
    # Annual count
    #--------------------------------
    if target=='track':
        a1bnd = np.arange(0,60+0.01,3) -0.5
    elif target=='point':
        #a1bnd = np.arange(0,0.15+0.001,0.02) -0.01
        a1bnd = np.arange(0,0.28+0.001,0.02) -0.01

    a1cnt = 0.5*(a1bnd[1:] + a1bnd[:-1])
    
    #-------------
    # Obs (annual count)
    #-------------
    a1num_obs = []
    for Year in ly_obs:
        if target =='track':
            stcid = set()
        elif target=='point':
            stcid = []

        for Mon in lmon:
            countdir = countbasedir + '/obs'
            util.mk_dir(countdir)
            obspath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)

            a1tcid_tmp = np.load(obspath)

            if target=='track':
                stcid = stcid.union(set(a1tcid_tmp))
            elif target=='point':
                stcid.extend(a1tcid_tmp)

        a1num_obs.append(len(stcid))

    if target=='point':
        nsteps = len(lmon)*30*4
        a1num_obs = np.array(a1num_obs) / float(nsteps)

    a1pdf_obs = np.histogram(a1num_obs, bins=a1bnd, density=True)[0]

    #a1bnd_tmp = a1bnd
    #a1pdf_obs_tmp = np.histogram(a1num_obs, bins=a1bnd_tmp)[0]
    #a1bnd_tmp = a1bnd_tmp[:-1]
    #print zip(a1bnd_tmp, a1pdf_obs_tmp)
    #-------------
    # Sim (annual count)
    #-------------
    d1pdf_sim = {}
    d2pdf_sim = {}
    for scen in lscen:
        a2pdf_sim = []

        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

        a1num_sim = []
        for ens in lens:
            for Year in ly_sim:
                if target =='track':
                    stcid = set()
                elif target=='point':
                    stcid = []
                for Mon in lmon:
                    countdir = countbasedir + '/%s/%s.%03d'%(slabel, scen, ens)
                    simpath = countdir + '/tcid.%04d.%02d.npy'%(Year,Mon)

                    try:
                        a1tcid_tmp = np.load(simpath)
                    except IOError:
                        a1tcid_tmp = np.array([])

                    if target=='track':
                        stcid = stcid.union(set(a1tcid_tmp))
                    elif target=='point':
                        stcid.extend(a1tcid_tmp)

                a1num_sim.append(len(stcid))

        if target=='point':
            nsteps = len(lmon)*30*4
            a1num_sim = np.array(a1num_sim) / float(nsteps)

        #-- Base histogram ----
        a1pdf_sim = np.histogram(a1num_sim, bins=a1bnd, density=True)[0]
        d1pdf_sim[scen] = a1pdf_sim

        #-- Range (bootstraping) ----------
        nset = 100
        nsample = len(a1num_sim)
        a2num_sim = np.array([np.random.choice(a1num_sim,nsample, replace=True) for i in range(nset)])

        a2pdf_sim = np.array([np.histogram(a2num_sim[i], bins=a1bnd, density=True)[0] for i in range(nset)])

        d2pdf_sim[scen] = a2pdf_sim
    #**************************************
    # Figure (Annual count histogram)
    #**************************************
    if figyear==True:
        fig = plt.figure(figsize=(7,4))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])

        a1pdf_his = d1pdf_sim['HPB']
        a1pdf_nat = d1pdf_sim['HPB_NAT']
        a1yerr_his_lw  =-np.percentile(d2pdf_sim['HPB'], 25, axis=0) + a1pdf_his
        a1yerr_his_up  = np.percentile(d2pdf_sim['HPB'], 75, axis=0) - a1pdf_his
        a1yerr_nat_lw  =-np.percentile(d2pdf_sim['HPB_NAT'], 25, axis=0) + a1pdf_nat
        a1yerr_nat_up  = np.percentile(d2pdf_sim['HPB_NAT'], 75, axis=0) - a1pdf_nat

        wbar = 0.25*(a1bnd[2] - a1bnd[1])
        a1x_obs = a1cnt - wbar
        a1x_nat = a1cnt
        a1x_his = a1cnt + wbar
        #ax.bar(a1x_obs, a1pdf_obs, width=wbar, color='0.3')
        #ax.bar(a1x_nat, a1pdf_nat, width=wbar, yerr=[a1yerr_nat_lw, a1yerr_nat_up], color='steelblue', ecolor='b')
        #ax.bar(a1x_his, a1pdf_his, width=wbar, yerr=[a1yerr_his_lw, a1yerr_his_up], color='tomato', ecolor='r')

        ax.plot(a1x_obs, a1pdf_obs, color='0.3')
        ax.plot(a1x_nat, a1pdf_nat, color='steelblue')
        ax.plot(a1x_his, a1pdf_his, color='tomato')



        xmax = {'track':30, 'point':a1bnd.max()}[target]
        plt.xlim([0,xmax])

        #plt.plot(a1bnd_tmp, a1pdf_obs_tmp)   # test

        stitle = target + '\n' + slabel + '\n' +sdomain
        plt.title(stitle)
        #plt.show()
        figpath = figdir + '/count.year.%s.%s.png'%(target, slabel)
        plt.savefig(figpath)
        print figpath



# %%

# %%
