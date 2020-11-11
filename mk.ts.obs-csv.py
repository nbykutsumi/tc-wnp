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

#ly_obs = range(2000,2010+1)
ly_obs = range(1980,2018+1)
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

for region in lregion:
    print region
    [[lllat,lllon],[urlat,urlon]] = dbbox[region]

    sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)
    
    outdir = '/home/utsumi/temp/bams2020/temp'
    util.mk_dir(outdir)
    
    #countbasedir = '/home/utsumi/temp/bams2020/count-%s/%s'%(target,sdomain)
    countbasedir = '/home/utsumi/temp/bams2020/count-point/%s'%(sdomain)



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

    print region, a2num_obs.shape
    sout = util.array2csv(a2num_obs)
    outpath = outdir + '/count.obs.%s.csv'%(region)
    f=open(outpath, 'w'); f.write(sout); f.close()
    print outpath
# %%
