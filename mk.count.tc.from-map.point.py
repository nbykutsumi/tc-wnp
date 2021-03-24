# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
#----------------------------------
import sys, os
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
import d4PDF
#--------------------------------------
figmon = True
#figmon = False
figyear= True
#figyear= False

#target='track'
target='point'

#ly_sim = range(2000,2010+1)
ly_sim = list(range(1980,2010+1))
#ly_obs = range(1989,1989+1)
#ly_obs = range(1980,2018+1)
#ly_obs = range(1980,2018+1)
#ly_obs = range(1980,2018+1)
ly_obs = list(range(1980,2010+1))
#ly_obs = range(2000,2010+1)
#lmon= range(1,12+1)
lmon= list(range(1,12+1))

#lregion = ['WNP', 'NEA','NNEA']
lregion = ['WNP']
dbbox = {
    #'WNP':[[20,120],[47,150]],
    #'NEA':[[30,120],[47,150]],
    'WNP':[[10,105],[47,150]],
    'NEA':[[20,115],[47,150]],
    'NNEA':[[25,120],[47,150]],
    }

dmax = {'WNP':1.0, 'NEA':0.5, 'NNEA':0.25}
dwbin= {'WNP':0.04, 'NEA':0.02, 'NNEA':0.01}
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
#lens    = list(range(1,19+1))
lens    = list(range(1,20+1))
#lens    = [1,2,3]
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

#**********************************
#----------------------------------
miss =-9999
miss_int= -9999


lthsst  = [27]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
#lthwcore= [-1,0,1]
lthdura = [36]
#lthwind = [10,13,15]
lthwind = [12]
#lthwdif = [-9999]
lthwdif = [-9999]



#lthsst  = [27,28]
#lexrvort = np.array([3])*1.0e-5
##ltcrvort = np.array([5])*1.0e-5
#ltcrvort = np.array([3])*1.0e-5
#lthwcore= [0]
##lthwcore= [-1,0,1]
#lthdura = [36]
##lthwind = [10,13,15]
#lthwind = [8,10,12,14]
##lthwdif = [-9999]
#lthwdif = [-9999,-30,-10]



lkey = [[thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif]
        for thsst   in lthsst
        for exrvort in lexrvort
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        for thdura in lthdura
        for thwind in lthwind
        for thwdif in lthwdif
        ]


for (thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif) in lkey:

    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)


    for region in lregion:
        print(region)
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
        for Year in ly_obs:
            a1num_tmp = []   # keep monthly data (12)
            for Mon in lmon:
                bstdir = outbaseDir + '/bst'

                #bstPath= bstdir + '/freq.tc.bst.%04d.%02d.pickle'%(Year,Mon)
                #with open(bstPath,'rb') as f:
                #    dbst = pickle.load(f)
                #a2count = dbst['a2count']
                #nstep   = dbst['nstep']
                a2count = np.load(bstdir + '/a2count.tc.bst.%04d.%02d.npy'%(Year,Mon))
                nstep   = np.load(bstdir + '/nstep.tc.bst.%04d.%02d.npy'%(Year,Mon))
   
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

                        print(region,scen,ens,Year,Mon)
                        #-- load --------
                        outDir = outbaseDir + '/%s/%s-%03d'%(slabel, scen, ens)
                        #objPath= outDir + '/freq.tc.obj.en-%03d.%04d.%02d.pickle'%(ens,Year,Mon)
                        #with open(objPath,'rb') as f:
                        #    dobj = pickle.load(f)
                        #a2count = dobj['a2count']
                        #nstep   = dobj['nstep']

                        a2count = np.load(outDir + '/a2count.tc.obj.%04d.%02d.npy'%(Year,Mon))
                        nstep   = np.load(outDir + '/nstep.tc.obj.%04d.%02d.npy'%(Year,Mon))


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
    
            fig = plt.figure(figsize=(7,4))
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
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
    
            stitle = target + '\n' + slabel + '\n' +sdomain
            plt.title(stitle)
            #plt.show()
            figpath = figdir + '/count.mon.%s.%s.%s.png'%(target, slabel, region)
            plt.savefig(figpath)
   
            print(figpath) 
    
        #********************************
        # Annual count
        #--------------------------------
        if target=='track':
            a1bnd = np.arange(0,60+0.01,3) -0.5
        elif target=='point':
            #a1bnd = np.arange(0,0.15+0.001,0.02) -0.01
            a1bnd = np.arange(0,dmax[region],dwbin[region]) -0.01
    
        a1cnt = 0.5*(a1bnd[1:] + a1bnd[:-1])
        
        #-------------
        # Obs (annual count)
        #-------------
        a1num_obs = a2num_obs.mean(axis=1)
        #a1pdf_obs = np.histogram(a1num_obs, bins=a1bnd, density=True)[0]
        a1pdf_obs = np.histogram(a1num_obs, bins=a1bnd)[0] / float(len(a1num_obs))
        #print a1num_obs 
        #-------------
        # Sim (annual count)
        #-------------
        d1pdf_sim = {}
        d2pdf_sim = {}

        d2pdf_sim3= {}
        for scen in lscen:
            a2pdf_sim = []
            a1num_sim = []

            a3num_sim = d3num_sim[scen]
            a2num_sim = a3num_sim.mean(axis=(2))
            a1num_sim = a2num_sim.flatten()
            #-- Base histogram ----
            #a1pdf_sim = np.histogram(a2num_sim, bins=a1bnd, density=True)[0]
            a1pdf_sim = np.histogram(a2num_sim, bins=a1bnd)[0] / float(len(a2num_sim.flatten()))
            d1pdf_sim[scen] = a1pdf_sim

            #-- Range (bootstraping) ----------
            nset = 100
            nsample = len(a1num_sim)
            a2num_sim = np.array([np.random.choice(a1num_sim,nsample, replace=True) for i in range(nset)])
    
            #a2pdf_sim = np.array([np.histogram(a2num_sim[i], bins=a1bnd, density=True)[0] for i in range(nset)])
            a2pdf_sim = np.array([np.histogram(a2num_sim[i], bins=a1bnd)[0] / float(len(a2num_sim[i])) for i in range(nset)])
    
            d2pdf_sim[scen] = a2pdf_sim




        #**************************************
        # Figure (Annual count histogram)
        #**************************************
        if figyear==True:
            fig = plt.figure(figsize=(7,4))
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    
            a1pdf_his = d1pdf_sim['HPB']
            #a1pdf_nat = d1pdf_sim['HPB_NAT']

            #a1yerr_his_lw  =-np.percentile(d2pdf_sim['HPB'], 25, axis=0) + a1pdf_his
            #a1yerr_his_up  = np.percentile(d2pdf_sim['HPB'], 75, axis=0) - a1pdf_his
            #a1yerr_nat_lw  =-np.percentile(d2pdf_sim['HPB_NAT'], 25, axis=0) + a1pdf_nat
            #a1yerr_nat_up  = np.percentile(d2pdf_sim['HPB_NAT'], 75, axis=0) - a1pdf_nat
   
            a1up_his = np.percentile(d2pdf_sim['HPB'], 95, axis=0)
            a1lw_his = np.percentile(d2pdf_sim['HPB'], 5, axis=0)

            #a1up_nat = np.percentile(d2pdf_sim['HPB_NAT'], 95, axis=0)
            #a1lw_nat = np.percentile(d2pdf_sim['HPB_NAT'], 5, axis=0)


            #-- Mask zero ---
            a1pdf_obs = ma.masked_equal(a1pdf_obs,0)
            #a1pdf_nat = ma.masked_equal(a1pdf_nat,0)
            a1pdf_his = ma.masked_equal(a1pdf_his,0)

            #a1up_nat  = ma.masked_equal(a1up_nat, 0)
            #a1lw_nat  = ma.masked_equal(a1lw_nat, 0)
            a1up_his  = ma.masked_equal(a1up_his, 0)
            a1lw_his  = ma.masked_equal(a1lw_his, 0)

            #----------------
 
            wbar = 0.25*(a1bnd[2] - a1bnd[1])
            a1x_obs = a1cnt - wbar
            #a1x_nat = a1cnt
            a1x_his = a1cnt + wbar
            ax.plot(a1x_obs, a1pdf_obs, linewidth=3,color='0.6')
            #ax.plot(a1x_nat, a1pdf_nat, linewidth=3,color='steelblue')
            ax.plot(a1x_his, a1pdf_his, linewidth=3,color='k')

            #ax.plot(a1x_nat, a1up_nat,'--', color='steelblue')
            #ax.plot(a1x_nat, a1lw_nat,'--', color='steelblue')
            ax.plot(a1x_his, a1up_his,'--', color='k')
            ax.plot(a1x_his, a1lw_his,'--', color='k')
    
            plt.yscale('log')    
            plt.ylim([0.001,1])   

 
            xmax = {'track':30, 'point':a1bnd.max()}[target]
            plt.xlim([0,xmax])
    
            #plt.plot(a1bnd_tmp, a1pdf_obs_tmp)   # test
    
            stitle = target + '\n' + slabel + '\n' +sdomain
            plt.title(stitle)
            #plt.show()
            figpath = figdir + '/count.year.%s.%s.%s.png'%(target, slabel, region)
            plt.savefig(figpath)
            print(figpath)



# %%

# %%
