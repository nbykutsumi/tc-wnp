# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
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
#ly_obs = range(1980,2018+1)
lly_obs = [range(1980,1999+1),range(2000,2010+1),range(2011,2018+1),range(1980,2018+1)]
#ly_obs = range(2000,2010+1)
#lmon= range(1,12+1)
lmon= range(1,12+1)

#lregion = ['WNP', 'NEA']
lregion = ['WNP']
dbbox = {
    'WNP':[[20,120],[47,150]],
    'NEA':[[30,120],[47,150]],
    }

dgrid    = 5.0
a1latcnt = np.arange(-90+ dgrid*0.5, 90, dgrid)
a1loncnt = np.arange(0+ dgrid*0.5, 360, dgrid)

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
lens    = [1,2]
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

#**********************************
#----------------------------------
miss =-9999
miss_int= -9999

a1lat = d4PDF.Lat()   # dlat ~ 0.5615674
a1lon = d4PDF.Lon()   # dlon = 0.5625

lthsst  = [27]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
#lthwcore= [-1,0,1]
lthdura = [36]
lthwind = [11]
#lthwind = [10,13,15]
#lthwdif = [-9999,-30,-20,-10]
lthwdif = [-9999]


lkey = [[thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif]
        for thsst   in lthsst
        for exrvort in lexrvort
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thdura in lthdura
        for thwind in lthwind
        for thwdif in lthwdif
        ]

for region in lregion:

    fig = plt.figure(figsize=(12,8))
    for iax,ly_obs in enumerate(lly_obs):

        ax = fig.add_subplot(2,2,iax+1)
        for (thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif) in lkey:

            exrvortout = exrvort*1.0e+5
            tcrvortout = tcrvort*1.0e+5

            #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
            slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)


            print region
            [[lllat,lllon],[urlat,urlon]] = dbbox[region]

            sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)

            outbaseDir = '/home/utsumi/temp/bams2020/map-freq'

            #countbasedir = '/home/utsumi/temp/bams2020/count-%s/%s'%(target,sdomain)
            countbasedir = '/home/utsumi/temp/bams2020/count-point/%s'%(sdomain)
            figdir  = '/home/utsumi/temp/bams2020/fig/temp'
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
            for Year in ly_obs:
                a1num_tmp = []   # keep monthly data (12)
                for Mon in lmon:
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
            #-------------
            # Sim (monthly count)
            #-------------
            d3num_sim = {}   
            for scen in lscen:
                a3num_sim = []  #  ENS x Y x M
                slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

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
            # Figure (mounthly count)
            #---------------------------

            if figmon is True:
                cm_ens = {'HPB':'0.5','HPB_NAT':'blue'}
                cm_ave = {'HPB':'tomato',  'HPB_NAT':'steelblue'}

                a1x = np.arange(1,12+1)

                a1num_obs_mean= a2num_obs.mean(axis=0)
                a1num_his_mean= d3num_sim['HPB'].mean(axis=(0,1))
                #a1num_nat_mean= d3num_sim['HPB_NAT'].mean(axis=(0,1))

                yerr_obs_lw =-np.percentile(a2num_obs, 25, axis=0) + a1num_obs_mean
                yerr_obs_up = np.percentile(a2num_obs, 75, axis=0) - a1num_obs_mean

                yerr_his_lw =-np.percentile(d3num_sim['HPB'], 25, axis=(0,1)) + a1num_his_mean

                yerr_his_up = np.percentile(d3num_sim['HPB'], 75, axis=(0,1)) - a1num_his_mean
                #yerr_nat_lw =-np.percentile(d3num_sim['HPB_NAT'], 25, axis=(0,1)) + a1num_nat_mean
                #yerr_nat_up = np.percentile(d3num_sim['HPB_NAT'], 75, axis=(0,1)) - a1num_nat_mean


                wbar = 0.25*(a1x[1] - a1x[0])
                a1x_nat = a1x - wbar
                a1x_his = a1x
                a1x_obs = a1x + wbar

                #ax.bar(a1x_nat, a1num_nat_mean, width=wbar, yerr=[yerr_nat_lw, yerr_nat_up], color='steelblue', ecolor='b')
                ax.bar(a1x_his, a1num_his_mean, width=wbar, yerr=[yerr_his_lw, yerr_his_up], color='k', ecolor='k')
                ax.bar(a1x_obs, a1num_obs_mean, width=wbar, yerr=[yerr_obs_lw, yerr_obs_up], color='0.3', ecolor='0.6')

                ax.set_ylim([0,0.8])
                iy,ey = ly_obs[0], ly_obs[-1]
                stitle = '%04d-%04d'%(iy,ey)
                plt.title(stitle)
    stitle = '%04d-%04d'%(iy,ey) + '\n' + slabel + '\n' +sdomain
    plt.suptitle(stitle)
    figpath = figdir + '/temp.count.mon.%s.%s.%s.png'%(target, slabel, region)
    plt.savefig(figpath)
    
    print figpath 


# %%
