# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
#----------------------------------
import sys, os, pickle
from   numpy import *
from   datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
import d4PDF
#--------------------------------------
#lens    = range(1,9+1)
#lens    = range(2,9+1)
lens    = range(1,4+1)

iYM = [2000,1]
eYM = [2000,12]

#iYM = [2000,1]
#eYM = [2010,12]
lYM = util.ret_lYM(iYM, eYM)
lYM = [(Year,Mon) for (Year,Mon) in lYM if Mon in [5,6,7,8,9,10]]

lscen   = ['HPB','HPB_NAT']
lvname = ['prec','wind']
#----------------------------------
a1latin = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonin = d4PDF.Lon()   # dlon = 0.5625

nyin    = len(a1latin)
nxin    = len(a1lonin)

miss =-9999
miss_int= -9999

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]
[[lllat,lllon],[urlat,urlon]] = [[20,120],[47,150]]
#[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]

for vname in lvname:
    d2hist_sim = {}
    for scen in lscen:
        d2hist_sim[scen] = []

        for ens in lens:
            a2hist_tmp = []
            for (Year,Mon) in lYM:
                #-- Read data -----
                datdir = '/home/utsumi/temp/bams2020/histogram/lat.%03d-%03d.lon.%04d-%04d/gen.%s.%03d'%(lllat,urlat,lllon,urlon, scen, ens)
                datpath = datdir + '/%04d.%02d.pickle'%(Year,Mon)
                with open(datpath, 'rb') as f:
                    ddat = pickle.load(f)

                a2hist_tmp.append(ddat[vname])
                a1bnd = ddat['bins-%s'%(vname)]

            a1hist_tmp = np.array(a2hist_tmp).sum(axis=0)
            d2hist_sim[scen].append(a1hist_tmp)
        d2hist_sim[scen] = np.array(d2hist_sim[scen])

    #** Make mean and range *****
    d1hist_mean = {}
    d1hist_up = {}
    d1hist_lw = {}

    for scen in lscen:
        d1hist_mean[scen] = d2hist_sim[scen].mean(axis=0)
        d1hist_up[scen] = np.percentile(d2hist_sim[scen], 90, axis=0)
        d1hist_lw[scen] = np.percentile(d2hist_sim[scen], 10, axis=0)


    #************************
    #---- Draw PDF ---------------
    fig = plt.figure(figsize=(6,5))
    ax  = fig.add_axes([0.1,0.1, 0.8, 0.8])

    #--- Simulation -----------
    for scen in lscen:
        cm_ens = {'HPB':'0.5','HPB_NAT':'blue'}[scen]
        cm_ave = {'HPB':'k',  'HPB_NAT':'royalblue'}[scen]

        a1cnt = 0.5*(a1bnd[1:] + a1bnd[:-1])
        plt.plot(a1cnt, d1hist_up[scen], '-', color=cm_ens, linewidth=0.8)
        plt.plot(a1cnt, d1hist_lw[scen], '-', color=cm_ens, linewidth=0.8)

        wbar= (a1bnd[2] - a1bnd[1])*0.4
        a1x = a1cnt + {'HPB':+0.5*wbar, 'HPB_NAT':-0.5*wbar}[scen]
        plt.bar(a1x, d1hist_mean[scen], color=cm_ave)

    stitle = 'Unconditional %s \n [[%d,%d],[%d,%d]]'%(vname, lllat,lllon,urlat,urlon)
    plt.title(stitle)
    if vname=='prec':
        plt.yscale('log')

    if vname=='prec':
        plt.xlim([0,60])
    elif vname=='wind':
        plt.xlim([0,60])

    plt.show()
    figdir = '/home/utsumi/temp/bams2020'
    figpath= figdir + '/pdf.gen.lat.%03d-%03d.lon.%04d-%04d.%s.png'%(lllat,lllon,urlat,urlon,vname)
    plt.savefig(figpath)
    print figpath


# %%
