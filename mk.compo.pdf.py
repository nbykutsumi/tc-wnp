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
#calcbst= True
calcbst= False
#figbst = True
figbst = False

calcobj= True
#calcobj= False
figobj = True
#figobj = False

figmean = True

#lens    = range(1,9+1)
lens    = [1]

ctype = 'TC'
#ctype = 'ALL'
iYM = [2000,8]
eYM = [2010,12]
lYM = util.ret_lYM(iYM,eYM)

cmbnd = [0,0.1, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50]
#cmbnd = None
dgridy = 9  # 9 x 0.5615674 ~ 5.05 degree radius
dgridx = 9  # 9 x 0.5625    ~ 5.06 degree radius
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
#lens    = [1]
res     = "320x640"
noleap  = False

detectName = 'wsd_d4pdf_20200428'
compbasedir= '/home/utsumi/temp/bams2020/composite'

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

    dprechist = {}
    dwindhist = {}        

    for scen in lscen:
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout, thwcore, thdura0)
        #slabel = 'wceach-%s-%s-%s-sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(ctype, expr, scen, thsst, exrvortout, tcrvortout, thwcore, thdura)


        for ens in lens:

            aprec = []
            awind = []
            for YM in lYM:
                Year,Mon = YM
                #--- Load file ---------
                compdir = compbasedir + '/%s/%s.%03d'%(slabel0,scen,ens)
                comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)
                with open(comppath,'rb') as f:
                    ddat = pickle.load(f)
                    print ddat
                    print comppath
                a3prec = ddat['prec']
                a3wind = ddat['wind']
                aprec.extend(a3prec.flatten().tolist())
                awind.extend(a3wind.flatten().tolist())
            plt.hist(awind)
            plt.savefig('/home/utsumi/temp/bams2020/temp.png')
            sys.exit() 

            #bndprec = np.arange(1,100+1,1)
            #cntprec = (bndprec[1:] + bndprec[:-1])*0.5
            #dprechist[scen,ens] = np.histogram(aprec, bins=bndprec, density=True)[0]
           
        #---- Draw PDF ---------------
        for ens in lens:
            plt.plot(cntprec, dprechist[scen,ens], '-',color='k')
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d'%(thsst0, exrvortout0, tcrvortout, thwcore, thdura0)
        figdir = '/home/utsumi/temp/bams2020'
        figpath= figdir + '/temp.pdf.prec.png'
        plt.savefig(figpath)
        print figpath
        
        
        
