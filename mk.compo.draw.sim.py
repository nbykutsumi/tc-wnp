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
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
#--------------------------------------
lens    = range(1,9+1)

iYM = [2000,1]
eYM = [2010,12]
lYM = util.ret_lYM(iYM,eYM)

#iYMobs=[1980,1]
#eYMobs=[2010,12]
iYMobs=[1980,1]
eYMobs=[1980,7]
lYMobs = util.ret_lYM(iYMobs, eYMobs)

a1bnd = (np.arange(0, 150+0.01, 1) -0.05) *1e-5
a1cent= 0.5* (a1bnd[:-1] + a1bnd[1:])
#varname = 'slp'
varname = 't500'
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
a3var_obs = []
for YM in lYMobs:
    Year,Mon = YM
    compdir  = compbasedir + '/obs' 
    comppath = compdir + '/%04d.%02d.pickle'%(Year,Mon)
    with open(comppath,'rb') as f:
        ddat = pickle.load(f)
    a3var_tmp = ddat[varname]
    if a3var_tmp.shape[0]==0: continue

    a3var_obs.append(a3var_tmp)

a3var_obs = np.concatenate(a3var_obs, axis=0)
a2var_obs = ma.masked_invalid(a3var_obs).mean(axis=0)
plt.imshow(a2var_obs)
plt.colorbar()
plt.show()
sys.exit()
# %%
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
    d1hist_mean_sim = {}
    for scen in lscen:
        a1hist_mean_sim = d2hist_sim[scen].sum(axis=0)
        a1hist_mean_sim = a1hist_mean_sim / float(a1hist_mean_sim.sum())
        d1hist_mean_sim[scen] = a1hist_mean_sim 

        #-- Normalize each ensemble ---
        d2hist_sim[scen] = d2hist_sim[scen]/d2hist_sim[scen].sum(axis=1).reshape(-1,1).astype('float32')

    #---- Draw PDF ---------------
    fig = plt.figure(figsize=(6,5))
    ax  = fig.add_axes([0.1,0.1, 0.8, 0.8])

    #--- Observation ----------
    plt.plot(a1cent, a1hist_obs, '-',color='r', linewidth=1.6)
    

    #--- Simulation -----------
    for scen in lscen: 
        cm_ens = {'HPB':'0.2','HPB_NAT':'royalblue'}[scen]
        cm_ave = {'HPB':'k',  'HPB_NAT':'blue'}[scen]
        
        for iens,ens in enumerate(lens):
            if scen=='HPB_NAT':continue
            plt.plot(a1cent, d2hist_sim[scen][iens], '-',color=cm_ens, linewidth=0.8)

        plt.plot(a1cent, d1hist_mean_sim[scen], '-', color=cm_ave, linewidth=1.6)
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout, thwcore, thdura0)
        plt.title(slabel)

        figdir = '/home/utsumi/temp/bams2020'
        figpath= figdir + '/pdf.vort.%s.png'%(slabel)
        plt.savefig(figpath)
        print figpath
            
        
        


# %%
