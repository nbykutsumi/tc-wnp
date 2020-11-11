# %%
import matplotlib
#matplotlib.use('Agg')
%matplotlib inline
#----------------------------------
import sys, os, pickle
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
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
figmon = True
#figmon = False
figyear= True
#figyear= False

#target='track'
target='point'

ly_sim = range(2000,2010+1)
#ly_obs = range(1989,1989+1)
#ly_obs = range(1980,2018+1)
#ly_obs = range(1980,2018+1)
ly_obs = range(2000,2010+1)
#ly_obs = range(1980,2018+1)
lmon= range(1,12+1)

lregion = ['WNP', 'NEA']
dbbox = {
    'WNP':[[20,120],[47,150]],
    'NEA':[[30,120],[47,150]],
    }

dgrid    = 5.0
a1latcnt = np.arange(-90+ dgrid*0.5, 90, dgrid)
a1loncnt = np.arange(0+ dgrid*0.5, 360, dgrid)

#-----------------
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
lens    = range(1,9+1)
#lens    = [1,2]
res     = "320x640"
noleap  = False

#**********************************
#----------------------------------
miss =-9999
miss_int= -9999

lthsst = [27]
exrvort= 3.0*1e-5
#ltcrvort= np.array([5, 6, 7])*1e-5
ltcrvort= np.array([5])*1e-5
lthwcore= [0]
#lthwcore= [3]
thdura = 36
lthwind = [15]
lthwdif = [-30]
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



for region in lregion:
    print region
    [[lllat,lllon],[urlat,urlon]] = dbbox[region]

    sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)
    
    outbaseDir = '/home/utsumi/temp/bams2020/map-freq'
    
    #countbasedir = '/home/utsumi/temp/bams2020/count-%s/%s'%(target,sdomain)
    countbasedir = '/home/utsumi/temp/bams2020/count-point/%s'%(sdomain)
    figdir  = '/home/utsumi/temp/bams2020/fig/%s'%(sdomain)
    util.mk_dir(figdir)



    a2loncnt, a2latcnt = np.meshgrid(a1loncnt, a1latcnt)
    a2masklat = ma.masked_outside(a2latcnt, lllat, urlat).mask
    a2masklon = ma.masked_outside(a2loncnt, lllon, urlon).mask
    a2maskregion = a2masklat + a2masklon

    #********************************
    # Monthly count
    #--------------------------------
    # Obs (monthly count)
    #-------------
    a2num_obs = []   # Y x M
    a1dtobs = []
    for Year in ly_obs:
        a1num_tmp = []   # keep monthly data (12)
        for Mon in lmon:
            a1dtobs.append(datetime(Year, Mon, 1))
            bstdir = outbaseDir + '/bst'
            bstPath= bstdir + '/freq.tc.bst.%04d.%02d.pickle'%(Year,Mon)
            with open(bstPath,'rb') as f:
                dbst = pickle.load(f)
   
            a2count = dbst['a2count']
            nstep   = dbst['nstep']
   
            avecount = ma.masked_where(a2maskregion, a2count).sum() / float(nstep)
            a1num_tmp.append(avecount)
    
        a2num_obs.append(a1num_tmp)
    
    a2num_obs = np.array(a2num_obs)

    ##----------------------------
    ## Figure (Monthly)
    ##----------------------------  
    #a1num_obs = a2num_obs.flatten()

    #fig = plt.figure(figsize=(12,3))
    #ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    #
    #ax.plot( a1dtobs, a1num_obs, color='0.5')
    #stitle = 'OBS %s'%(region)
    #ax.grid()
    #plt.title(stitle)

    ##----------------------------
    ## Figure (Annual)
    ##----------------------------  
    #a1num_obs = a2num_obs.mean(axis=1)
    #a1x       = [datetime(Y,1,1) for Y in ly_obs]

    #fig = plt.figure(figsize=(12,3))
    #ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    #
    #ax.plot( a1x, a1num_obs, color='0.5')
    #stitle = 'OBS %s'%(region)
    #plt.title(stitle)
    #ax.xaxis.set_major_locator(mdates.YearLocator(5))
    #ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ##ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    #ax.grid()
    #ax.grid(which='minor', linestyle='--')


    ## %%

    #-------------
    # Sim (load data)
    #-------------

    for (thsst,tcrvort,thwcore,thwind,thwdif) in lkey:
    
        exrvortout = exrvort*1.0e+5
        tcrvortout = tcrvort*1.0e+5
    
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)


        d3num_sim = {}   
        for scen in lscen:
            a3num_sim = []  #  ENS x Y x M
            slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    
            for ens in lens:
                for Year in ly_sim:

                    a1num_tmp = []
                    for Mon in lmon:

                        #-- load --------
                        outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                        objPath= outDir + '/freq.tc.obj.en-%03d.%04d.%02d.pickle'%(ens,Year,Mon)
                        with open(objPath,'rb') as f:
                            dobj = pickle.load(f)
                        

                        a2count = dobj['a2count']
                        nstep   = dobj['nstep']
                        avecount = ma.masked_where(a2maskregion, a2count).sum() / float(nstep)
                        a3num_sim.append(avecount)

            a3num_sim = np.array(a3num_sim).reshape(len(lens),len(ly_sim),-1)
    
            d3num_sim[scen] = a3num_sim

        #---------------------------
        # Figure (annual count)
        #---------------------------
        a1num_obs = a2num_obs.mean(axis=1)
        a1x_obs   = [datetime(Y,1,1) for Y in ly_obs]

        a2num_his = d3num_sim['HPB'].mean(axis=2)
        a1num_his = d3num_sim['HPB'].mean(axis=(0,2))
        a1x_sim   = [datetime(Y,1,1) for Y in ly_sim]

        fig = plt.figure(figsize=(9,3))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        
        ax.plot( a1x_obs, a1num_obs, color='0.5', linewidth=2)
        ax.plot( a1x_sim, a1num_his, color='k', linewidth=2)
        for i in range(len(lens)):
            ax.plot( a1x_sim, a2num_his[i], color='0.7', linewidth=0.5)


        stitle = '%s %s'%(region, slabel)
        plt.title(stitle)
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        ax.xaxis.set_minor_locator(mdates.YearLocator(1))
        #ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.grid()
        ax.grid(which='minor', linestyle='--')

        #---------------------------
        # Figure (monthly count)
        #---------------------------
        a1num_obs = a2num_obs.flatten()
        a1x_obs   = a1dtobs

        a2num_his = d3num_sim['HPB'].reshape(len(lens), -1)
        a1num_his = d3num_sim['HPB'].mean(axis=0).flatten()
        a1x_sim   = [datetime(Y,M,1) for Y in ly_sim for M in range(1,12+1)]

        fig = plt.figure(figsize=(18,6))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        
        ax.plot( a1x_obs, a1num_obs, color='0.5', linewidth=2)
        ax.plot( a1x_sim, a1num_his, color='k', linewidth=2)
        for i in range(len(lens)):
            ax.plot( a1x_sim, a2num_his[i], color='0.7', linewidth=0.5)


        stitle = '%s %s'%(region, slabel)
        plt.title(stitle)
        ax.xaxis.set_major_locator(mdates.YearLocator(5))
        ax.xaxis.set_minor_locator(mdates.YearLocator(1))
        #ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax.grid()
        ax.grid(which='minor', linestyle='--')

       

      


# %%
