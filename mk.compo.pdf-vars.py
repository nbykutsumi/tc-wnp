# %%
import matplotlib
#matplotlib.use('Agg')
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
#mode = "validation"
mode = "difference"

prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#scen    = 'HPB' # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT']
if mode=="validation":
    lscen   = ['HPB']
    iYM = [1980,1]
    eYM = [2002,12]

elif mode=="difference":
    lscen   = ['HPB',"HPB_NAT"]
    iYM = [1990,1]
    eYM = [2010,12]



else:
    print("check mode",mode)
    sys.exit()
lens    = range(1,9+1)
res     = "320x640"
noleap  = False

lYM = util.ret_lYM(iYM,eYM)

iYMobs=[1980,1]
eYMobs=[2002,12]
lYMobs = util.ret_lYM(iYMobs, eYMobs)

#lvname = ['prec','wind','slp']
#lvname = ['prec','slp']
#lvname = ['prec-06']
lvname = ['prec000']
#lvname = ['wind']
lregion= ["WNP"]
dbbox ={
    'WNP':[[10,105],[47,150]],
    'NEA':[[20,115],[47,150]],
    'NNEA':[[25,120],[47,150]],
}


detectName = 'wsd_d4pdf_20200428'
compbasedir= '/home/utsumi/temp/bams2020/composite'
nradorg = 9  # composit: (2*nradorg+1)*(2*nradorg+1)

for region in lregion:
    [[lllat,lllon],[urlat,urlon]] = dbbox[region]

    for vname in lvname:
        if (len(vname)>=4)&(vname[:4]=='prec'):
            #a1bnd = np.arange(0, 80+0.01, 1)
            a1bnd = np.arange(0, 60.01, 2.5)
        elif vname=='wind':
            a1bnd = np.arange(0,50,2)
        elif vname=='slp':
            a1bnd = np.arange(930,1010+0.01,2)

        a1cent= 0.5* (a1bnd[:-1] + a1bnd[1:])
        #**********************************
        # Obs
        #----------------------------------
        a1var_obs = []
        for YM in lYMobs:
            Year,Mon = YM
            compdir  = compbasedir + '/obs-era5-tc/%04d'%(Year)
            a3var = np.load(compdir + "/%s.%04d.%02d.npy"%(vname,Year,Mon))
            a1lat = np.load(compdir + "/lat.%04d.%02d.npy"%(Year,Mon))
            a1lon = np.load(compdir + "/lon.%04d.%02d.npy"%(Year,Mon))

            if a3var.shape[0]==0: continue


            #-- Screening region ----
            a1flaglat = ma.masked_inside(a1lat, lllat, urlat).mask
            a1flaglon = ma.masked_inside(a1lon, lllon, urlon).mask

            a1flag = a1flaglat * a1flaglon
            #------------------------
            a1var_tmp = a3var[a1flag,:,:].flatten()
            #a1var_tmp = a3var[a1flag,:,:].max(axis=(1,2))
            a1var_obs.extend(a1var_tmp.tolist())

        a1var_obs  = np.array(a1var_obs)

        #-- Unit conversion -----
        if vname=='slp':
            a1var_obs = a1var_obs * 0.01

        if (len(vname)>=4)&(vname[:4]=='prec'):
            a1var_obs = a1var_obs *1000.   # m/h --> mm/h

        a1hist_obs, _ = np.histogram(a1var_obs, bins=a1bnd)
        a1hist_obs = a1hist_obs / float(a1hist_obs.sum())
        #----------------------------------
        # Input parameters
        #----------------------------------
        miss =-9999
        miss_int= -9999
        thsst = -9999
        exrvortout = 3
        tcrvortout = 3
        thwcore  = 0
        thwind   = -9999
        thwdif   = -9999  # not used
        thdura   = 36
        slabel0 = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)


        #--- Output parapeters ---
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

        slabel_out = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
        #--------------------------

        d2_sim = {}
        d2hist_sim = {}
        for scen in lscen:
            d2hist_sim[scen] = []

            for ens in lens:
                a1var_sim = []
                for YM in lYM:                           
                    Year,Mon = YM
                    #--- Load file ---------
                    compdir = compbasedir + '/%s/%s.%03d/%04d'%(slabel0,scen,ens,Year)
                    a3var = np.load(compdir + "/%s.%04d.%02d.npy"%(vname,Year,Mon))
                    a1vort= np.load(compdir + "/%s.%04d.%02d.npy"%("vort",Year,Mon))
                    a1tsum= np.load(compdir + "/%s.%04d.%02d.npy"%("tsum",Year,Mon))
                    a1wlw = np.load(compdir + "/%s.%04d.%02d.npy"%("wlw",Year,Mon))
                    a1lat = np.load(compdir + "/%s.%04d.%02d.npy"%("lat",Year,Mon))
                    a1lon = np.load(compdir + "/%s.%04d.%02d.npy"%("lon",Year,Mon))

                    if a3var.shape[0] ==0: continue

                    #-- Screen TC -----------
                    a1flagvort = ma.masked_greater_equal(a1vort, tcrvort).mask
                    a1flagwc   = ma.masked_greater_equal(a1tsum, thwcore).mask
                    a1flagwind = ma.masked_greater_equal(a1wlw, thwind).mask
                    #-- Region mask ---------
                    a1flaglat = ma.masked_inside(a1lat, lllat, urlat).mask
                    a1flaglon = ma.masked_inside(a1lon, lllon, urlon).mask

                    #------------------------
                    a1flag = a1flagvort * a1flagwc * a1flagwind * a1flaglat* a1flaglon

                    if a1flag is np.bool_(False):
                        continue

                    if vname=='slp':
                        a1var_tmp = a3var[a1flag,nradorg,nradorg]*0.01
        
                    elif (len(vname)>=4)&(vname[:4]=='prec'):
                        a1var_tmp = a3var[a1flag,nradorg,nradorg]*60*60

                    else:
                        a1var_tmp = a3var[a1flag,:,:].flatten()
                        #a1var_tmp = a3var[a1flag,:,:].max(axis=(1,2))

                    a1var_sim.extend(a1var_tmp.tolist())

                a1var_sim = np.array(a1var_sim)
                a1hist_sim, _ = np.histogram(a1var_sim, bins=a1bnd)

                d2hist_sim[scen].append(a1hist_sim)

            d2hist_sim[scen] = np.array(d2hist_sim[scen])

        #---- Make mean and normalize simulation --
        d1hist_mean_sim = {}
        for scen in lscen:
            a1hist_mean_sim = d2hist_sim[scen].sum(axis=0)
            d1hist_mean_sim[scen] = a1hist_mean_sim / float(a1hist_mean_sim.sum())
            #a1hist_mean_sim = a1hist_mean_sim / float(a1hist_mean_sim.sum())
            #d1hist_mean_sim[scen] = a1hist_mean_sim / float(len(lens))

            #-- Normalize each ensemble ---
            #d2hist_sim[scen] = d2hist_sim[scen]/d2hist_sim[scen].sum(axis=1).reshape(-1,1).astype('float32')

        #---- Range ----
        d1hist_up = {}
        d1hist_lw = {}
        d1err_up  = {}
        d1err_lw  = {}
        for scen in lscen:
            d1hist_up[scen] = np.percentile(d2hist_sim[scen], 75, axis=0)
            d1hist_lw[scen] = np.percentile(d2hist_sim[scen], 25, axis=0)
            d1err_up[scen]  = d1hist_up[scen] - d1hist_mean_sim[scen]
            d1err_lw[scen]  = -d1hist_lw[scen] + d1hist_mean_sim[scen]
        #---- Draw PDF ---------------
        fig = plt.figure(figsize=(10,4))
        ax  = fig.add_axes([0.1,0.1, 0.8, 0.8])

        #--- Observation ----------
        plt.plot(a1cent, a1hist_obs, '-',color='gray', linewidth=1.4)


        #--- Simulation -----------
        for scen in lscen: 
            cm_ens = {'HPB':'0.5','HPB_NAT':'blue'}[scen]
            cm_ave = {'HPB':'black',  'HPB_NAT':'steelblue'}[scen]

            #plt.plot(a1cent, d1hist_up[scen], '-', color=cm_ens, linewidth=0.8)
            #plt.plot(a1cent, d1hist_lw[scen], '-', color=cm_ens, linewidth=0.8)

            plt.plot(a1cent, d1hist_mean_sim[scen], '-', color=cm_ave, linewidth=2.5)
            #wbar= (a1bnd[2] - a1bnd[1])*0.4
            #a1x = a1cent + {'HPB':+0.5*wbar, 'HPB_NAT':-0.5*wbar}[scen]
            #plt.bar(a1x, d1hist_mean_sim[scen], color=cm_ave)
            #plt.bar(a1x, d1hist_mean_sim[scen], width=wbar, color=cm_ave, yerr=[d1err_lw[scen], d1err_up[scen]], ecolor=cm_ens)

        stitle = vname + ' ' + slabel_out
        stitle = stitle + '\n' + '[[%d,%d],[%d,%d]]'%(lllat,lllon,urlat,urlon)
        stitle = stitle + " %04d-%04d"%(iYM[0],eYM[0])
        plt.title(stitle)
        if (len(vname)>=4)&(vname[:4]=='prec'):
            plt.yscale('log')

        #plt.xlim([30,50])
        #plt.ylim([1e-6,1e-1])
        plt.show()            
        figdir = '/home/utsumi/temp/bams2020'
        figpath= figdir + '/pdf.%s.%s.%04d-%04d.png'%(vname,slabel_out,iYM[0],eYM[0])
        plt.savefig(figpath)
        print(figpath)




# %%
