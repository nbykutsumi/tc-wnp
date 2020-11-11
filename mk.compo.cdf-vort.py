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
lens    = range(1,9+1)

iYM = [2000,1]
eYM = [2010,12]
#eYM = [2001,12]
lYM = util.ret_lYM(iYM,eYM)

iYMobs=[1980,1]
eYMobs=[2010,12]
lYMobs = util.ret_lYM(iYMobs, eYMobs)

a1bnd = (np.arange(0, 500+0.01, 1) -0.05) *1e-5
a1cent= 0.5* (a1bnd[:-1] + a1bnd[1:])

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
#lscen   = ['HPB_NAT']
#lens    = range(1,9+1)
#lens    = [1]
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
compbasedir= '/home/utsumi/temp/bams2020/composite'

#**********************************
# Obs
#----------------------------------
a1vort_obs = []
for YM in lYMobs:
    Year,Mon = YM
    compdir  = compbasedir + '/obs' 
    comppath = compdir + '/%04d.%02d.pickle'%(Year,Mon)
    with open(comppath,'rb') as f:
        ddat = pickle.load(f)
    a1vort_tmp = ddat['vort']
    a1vort_obs.extend(a1vort_tmp.tolist())

a1vort_obs = np.array(a1vort_obs)
a1hist_obs, _ = np.histogram(a1vort_obs, bins=a1bnd)

a1cdf_obs = np.cumsum(a1hist_obs)
a1cdf_obs = a1cdf_obs / float(a1hist_obs.sum())
#----------------------------------
miss =-9999
miss_int= -9999
thsst0 = 27
exrvortout0 = 3
tcrvortout0 = 3
thwcore0 = 0
thdura0  = 48

slabel0 = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout0, thwcore0, thdura0)

ltcrvort = np.array([4.8])*1.0e-5
lthwcore = np.array([0.3])

lkey = [[tcrvort,thwcore]
        for tcrvort in ltcrvort
        for thwcore in lthwcore
        ]
print 'AAAAAAAAAAAAA'
for (tcrvort,thwcore) in lkey:
    tcrvortout = tcrvort*1.0e+5


    d2vort_sim = {}
    d2hist_sim = {}
    for scen in lscen:
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout, thwcore, thdura0)
        #slabel = 'wceach-%s-%s-%s-sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(ctype, expr, scen, thsst, exrvortout, tcrvortout, thwcore, thdura)
        d2vort_sim[scen] = []
        d2hist_sim[scen] = []
        for ens in lens:
            a1vort_sim = []
            for YM in lYM:
                Year,Mon = YM
                #--- Load file ---------
                compdir = compbasedir + '/%s/%s.%03d'%(slabel0,scen,ens)
                comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)
                with open(comppath,'rb') as f:
                    ddat = pickle.load(f)
                a1vort_tmp = ddat['vort']
                a1histsim, a1bnd = np.histogram(a1vort_tmp, bins=a1bnd)

                a1vort_sim.extend(a1vort_tmp.tolist())

            a1vort_sim = np.array(a1vort_sim)
            a1hist_sim, _ = np.histogram(a1vort_sim, bins=a1bnd)

            d2vort_sim[scen].append(a1vort_sim)
            d2hist_sim[scen].append(a1hist_sim)

        d2vort_sim[scen] = np.array(d2vort_sim[scen])
        d2hist_sim[scen] = np.array(d2hist_sim[scen])

    #---- Make mean and normalize simulation --
    d1cdf_all_sim = {}
    for scen in lscen:
        a1hist_all_sim = d2hist_sim[scen].sum(axis=0)

        a1cdf_all_sim = np.cumsum(a1hist_all_sim)
        a1cdf_all_sim = a1cdf_all_sim / float(a1hist_all_sim.sum())
        d1cdf_all_sim[scen] = a1cdf_all_sim
    #---- Draw CDF ---------------
    fig = plt.figure(figsize=(6,5))
    ax  = fig.add_axes([0.1,0.1, 0.8, 0.8])

    #--- Observation ----------
    print len(a1cent), d2hist_sim['HPB'].shape
    plt.plot(a1cent, a1cdf_obs, '-',color='r', linewidth=1.6)
    

    #--- Simulation -----------
    for scen in lscen: 
        cm_ave = {'HPB':'k', 'HPB_NAT':'b'}[scen]
        print len(a1cent)
        plt.plot(a1cent, d1cdf_all_sim[scen], '-', color=cm_ave, linewidth=1.6)
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout, thwcore, thdura0)
        plt.title(slabel)

    plt.xlim([0,0.002])
    figdir = '/home/utsumi/temp/bams2020'
    figpath= figdir + '/cdf.vort.%s.png'%(slabel)
    plt.savefig(figpath)
    print figpath
            
        
        


# %%
