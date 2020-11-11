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
#lens    = range(1,9+1)
#lens    = range(2,9+1)
lens    = [1]

iYM = [2000,1]
eYM = [2010,12]

#iYM = [2000,9]
#eYM = [2000,9]

lYM = util.ret_lYM(iYM,eYM)

#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
lscen   = ['HPB_NAT']
#lens    = range(1,9+1)
#lens    = [1]

figdir  = '/home/utsumi/temp/bams2020'
compbasedir= '/home/utsumi/temp/bams2020/composite'
#----------------------------------
miss =-9999
miss_int= -9999

#************************************
# d4PDF (Objective detection)
#************************************
thsst    = 27
exrvort  = 3*1.e-5
tcrvort  = 3*1.e-5
thwcore  = 0
thdura   = 48
exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

for scen in lscen:
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thdura)

    for ens in lens:
        print 'ens=',ens
        for YM in lYM:
            Year,Mon = YM
            lDTime = util.ret_lDTime_fromYM([Year,Mon],[Year,Mon], timedelta(hours=6), hour0=0)

            #--- Load composite pickle file ---------
            compdir = compbasedir + '/%s/%s.%03d'%(slabel,scen,ens)
            comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)
            with open(comppath,'rb') as f:
                 
                dat = pickle.load(dout, f)
             
            a3wind  = ddat['wind']
            a1idate = ddat['idate']
            a1ipos  = ddat['ipos']
            

# %%
