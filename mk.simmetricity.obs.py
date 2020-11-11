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
from   scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
#--------------------------------------
#iYM = [1980,1]
#eYM = [1980,1]

iYM = [2010,1]
eYM = [2013,12]

iYM_tc = [1980,1]
eYM_tc = [2018,12]

lYM    = util.ret_lYM(iYM,eYM)
lYM_tc = util.ret_lYM(iYM_tc,eYM_tc)

calc   = False
#calc   = True
lctype = ['tc','exc']
#lctype = ['exc']
#lctype = ['tc']

lvarname = ['t850']
#lvarname = ['t500']
#lvarname = ['prec']
#lvarname = ['slp']
#lvarname = ['wind']
#lvarname = ['slp', 'prec', 'wind850', 'wind500','wind', 't500', 't850']

nrad = 4
#-----------------
compbasedir= '/home/utsumi/temp/bams2020/composite'

def symmetricity(a):
    cc1 = np.ma.corrcoef(a.flatten(),np.rot90(a,1).flatten())[0,1]
    cc2 = np.ma.corrcoef(a.flatten(),np.rot90(a,2).flatten())[0,1]
    #cc1 = np.corrcoef(a.flatten(),a[::-1,:].flatten())[0,1]
    #cc2 = np.corrcoef(a.flatten(),a[:,::-1].flatten())[0,1]


    cc0 = (cc1 + cc2)*0.5
    return cc0,cc1,cc2
#**********************************
# Obs
#----------------------------------
for varname in lvarname:
    for ctype in lctype:
        if calc is False: continue

        if ctype=='exc':
            lYM_tmp = lYM
        elif ctype=='tc':
            lYM_tmp = lYM_tc
        else:
            print 'check', ctype
            sys.exit()

        a3var_obs = []
        a3slp_obs = []  # test
        a1lat_obs = []
        a1lon_obs = []
        lcc0      = []
        lcc1      = []
        lcc2      = []
        llat      = []
        llon      = []
        for YM in lYM_tmp:
            Year,Mon = YM
            compdir  = compbasedir + '/obs-%s'%(ctype)
            comppath = compdir + '/%04d.%02d.pickle'%(Year,Mon)
            with open(comppath,'rb') as f:
                ddat = pickle.load(f)
            a3var_tmp = ddat[varname]
            a3slp_tmp = ddat['slp']
            a1lat_tmp = ddat['lat']
            a1lon_tmp = ddat['lon']

            if a3var_tmp.shape[0]==0: continue

            a3var_tmp = a3var_tmp[:,9-nrad:9+nrad+1,9-nrad:9+nrad+1]
            a3slp_tmp = a3slp_tmp[:,9-nrad:9+nrad+1,9-nrad:9+nrad+1]

            a3var_obs.append(a3var_tmp)
            a3slp_obs.append(a3slp_tmp)  # tesmp
            a1lat_obs.append(a1lat_tmp)
            a1lon_obs.append(a1lon_tmp)

        a3var_obs = np.concatenate(a3var_obs, axis=0)
        a3slp_obs = np.concatenate(a3slp_obs, axis=0)
        a1lat_obs = np.concatenate(a1lat_obs)
        a1lon_obs = np.concatenate(a1lon_obs)

        a1mean_obs= ma.masked_invalid(a3var_obs).mean(axis=(1,2))
        a2var_obs = ma.masked_invalid(a3var_obs - a1mean_obs.reshape(-1,1,1)).mean(axis=0)

        #-- test -------
        ymin,xmin = np.unravel_index(np.argmin(a2var_obs), a2var_obs.shape)

        print 'shape_3d=',a3var_obs.shape
        print 'shape_2d=',a2var_obs.shape
        print 'ymin, xmin=', ymin, xmin
        #** Figure *******
        fig = plt.figure(figsize=(4,4))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])
        plt.imshow(a2var_obs, origin='lower')

        plt.plot(xmin,ymin,'o','r')

        plt.colorbar()
        stitle = '%s %s'%(ctype, varname)
        plt.title(stitle)
        plt.show()
        figdir = '/home/utsumi/temp/bams2020/composite/fig'
        figpath= figdir + '/composite.%s.%s.png'%(ctype,varname)
        print figpath 
#        sys.exit()

        #--- For each ------------
        for i,a in enumerate(a3var_obs):
            a2slp = a3slp_obs[i]
            ymin,xmin = np.unravel_index(np.argmin(a2slp), a2slp.shape)
            if (ymin !=nrad)or(xmin !=nrad): continue

            cc0,cc1,cc2 = symmetricity(a)
            lcc0.append(cc0)
            lcc1.append(cc1)
            lcc2.append(cc2)
            llat.append(a1lat_obs[i])
            llon.append(a1lon_obs[i])
            #print cc1,cc2


            ##-- test figure ------
            #print ddat['time'][i], ddat['lat'][i],ddat['lon'][i]
            #print '%.4f,%.4f,%.4f'%(lcc)
            #plt.imshow(a, origin='lower')
            ##plt.plot(xmin,ymin,'o','r')
            #plt.plot(nrad,nrad,'o',color='r')
            #plt.plot(xmin,ymin,'x',color='r',markersize=10)
            #plt.colorbar()
            #plt.title('%s %s %d'%(ctype, varname, i))
            #plt.show()

        #-- save data ---
        lcc0 = np.array(lcc0)
        lcc1 = np.array(lcc1)
        lcc2 = np.array(lcc2)
        llat = np.array(llat)
        llon = np.array(llon)
        outdir = '/home/utsumi/temp/bams2020/composite/npy'
        outpath0 = outdir + '/symmetricity.cc0.%s.%s.npy'%(ctype,varname)      
        outpath1 = outdir + '/symmetricity.cc1.%s.%s.npy'%(ctype,varname)      
        outpath2 = outdir + '/symmetricity.cc2.%s.%s.npy'%(ctype,varname)      
        latpath  = outdir + '/lat.%s.npy'%(ctype) 
        lonpath  = outdir + '/lon.%s.npy'%(ctype) 

        np.save(outpath0, lcc0)
        np.save(outpath1, lcc1)
        np.save(outpath2, lcc2)
        np.save(latpath,  llat)
        np.save(lonpath,  llon)
        print outpath0

    #*********************
    # PDF
    #*********************
    for i in range(3):
        outdir = '/home/utsumi/temp/bams2020/composite/npy'
        a1tc  = np.load(outdir + '/symmetricity.cc%d.tc.%s.npy'%(i,varname) )
        a1ex  = np.load(outdir + '/symmetricity.cc%d.exc.%s.npy'%(i,varname) )

        a1bin = np.arange(-1,1+0.001, 0.05)
        a1cnt = 0.5*(a1bin[1:]+a1bin[:-1])
        a1hist_tc,_ = np.histogram(a1tc, bins=a1bin, density=True)
        a1hist_ex,_ = np.histogram(a1ex, bins=a1bin, density=True)

        plt.figure(figsize=(3,3))
        plt.plot(a1cnt, a1hist_tc,'-',color='k')
        plt.plot(a1cnt, a1hist_ex,'--',color='k')
        plt.title(i)
        plt.show()

    #*********************
    # CDF
    #*********************
    for i in range(3):
        outdir = '/home/utsumi/temp/bams2020/composite/npy'
        a1tc  = np.load(outdir + '/symmetricity.cc%d.tc.%s.npy'%(i,varname) )
        a1ex  = np.load(outdir + '/symmetricity.cc%d.exc.%s.npy'%(i,varname) )

        a1bin = np.arange(-1,1+0.001, 0.05)
        a1cnt = 0.5*(a1bin[1:]+a1bin[:-1])
        a1hist_tc,_ = np.histogram(a1tc, bins=a1bin)
        a1hist_ex,_ = np.histogram(a1ex, bins=a1bin)

        a1cdf_tc = np.cumsum(a1hist_tc) / float(len(a1tc))
        a1cdf_ex = np.cumsum(a1hist_ex) / float(len(a1ex))

        plt.figure(figsize=(3,3))
        plt.plot([0,0],[0,1],'-')
        plt.plot(a1cnt, a1cdf_tc,'-', color='k')
        plt.plot(a1cnt, a1cdf_ex,'--',color='k')
        plt.title(i)
        plt.show()
    #*********************
    # Latitude & simmetricity
    #*********************
    for i in range(3):
        outdir = '/home/utsumi/temp/bams2020/composite/npy'
        a1tc  = np.load(outdir + '/symmetricity.cc%d.tc.%s.npy'%(i,varname) )
        a1ex  = np.load(outdir + '/symmetricity.cc%d.exc.%s.npy'%(i,varname) )

        a1lattc= np.load(outdir + '/lat.tc.npy')
        a1latex= np.load(outdir + '/lat.exc.npy')

        a1bin = np.arange(20,45,2.5)
        a1cnt = 0.5*(a1bin[1:] + a1bin[:-1])

        a1avetc, _, _ = stats.binned_statistic(x=a1lattc, values=a1tc, statistic='mean', bins=a1bin)
        a1aveex, _, _ = stats.binned_statistic(x=a1latex, values=a1ex, statistic='mean', bins=a1bin)

        a1numtc, _, _ = stats.binned_statistic(x=a1lattc, values=a1tc, statistic='count', bins=a1bin)
        a1numex, _, _ = stats.binned_statistic(x=a1latex, values=a1ex, statistic='count', bins=a1bin)

        plt.figure(figsize=(3,3))
        plt.title('simmetricity'+'%d'%i)
        plt.plot(a1cnt, a1avetc, '-', color='k')
        plt.plot(a1cnt, a1aveex, '--', color='k')
        plt.show()

        plt.figure(figsize=(3,3))
        plt.title('simmetricity'+'%d'%i)
        plt.plot(a1latex, a1ex,'.',color='b')
        plt.plot(a1lattc, a1tc,'.',color='k')
        plt.xlim([20,45])
        plt.show()

        plt.figure(figsize=(3,3))
        plt.title('count'+'%d'%i)
        plt.plot(a1cnt, a1numtc, '-', color='k')
        plt.plot(a1cnt, a1numex, '--', color='k')
        plt.show()

    #*********************
    # simmetricity 1 & 2
    plt.figure(figsize=(3,3))
    for ctype in ['tc','exc']:
        outdir = '/home/utsumi/temp/bams2020/composite/npy'
        a1var1  = np.load(outdir + '/symmetricity.cc1.%s.%s.npy'%(ctype,varname) )
        a1var2  = np.load(outdir + '/symmetricity.cc2.%s.%s.npy'%(ctype,varname) )

        if ctype=='tc': mycm='k'
        elif ctype=='exc':mycm='b'
        plt.plot(a1var1, a1var2,'.',color=mycm)
    plt.title('simmetricity 1 vs 2'
    plt.show()
   #*********************


# %%
