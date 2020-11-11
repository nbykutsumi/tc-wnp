# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
#----------------------------------
import sys, os
from   numpy import *
from   datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import util
from bisect import bisect_left
from detect_fsub import *
import d4PDF
#--------------------------------------
lens    = range(1,9+1)
#lens    = range(2,9+1)
#lens    = [1]

lvname  = ['prec','wind']
#lvname  = ['wind']
lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB']
lYear = range(2000,2010+1)
#lYear = range(2000,2001+1)
lMon  = range(1,12+1)
d4pdfdir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
d4atm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir=d4pdfdir)

a1latin = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonin = d4PDF.Lon()   # dlon = 0.5625

nyin    = len(a1latin)
nxin    = len(a1lonin)

miss =-9999
miss_int= -9999

lkeys = [[vname,scen,ens] for vname in lvname for scen in lscen for ens in lens]
for (vname,scen,ens) in lkeys:

    #if (vname=='prec')&(scen=='HPB'):continue  # test

    a3max = []
    for Year in lYear:
        a3tmp = []
        for Mon in lMon:
            if vname=='prec':
                a3in = d4sfc.load_6hr_mon('PRECIPI',scen,ens,Year,Mon)
            elif vname=='wind':
                a3u  = ma.masked_equal(d4sfc.load_6hr_mon('UAOPN',scen,ens,Year,Mon), -9999.)
                a3v  = ma.masked_equal(d4sfc.load_6hr_mon('VAOPN',scen,ens,Year,Mon), -9999.)
                a3in = np.sqrt( np.square(a3u) + np.square(a3v) )
            else:
                print 'check vname',vname
                sys.exit()

            #-- find max --
            a2tmp = a3in.max(axis=0)
            print vname,scen,ens,Year,Mon
            a3tmp.append(a2tmp)

        a3tmp = np.array(a3tmp)
        a2max = a3tmp.max(axis=0)
        a3max.append(a2max)
    a3max = np.array(a3max)
    
    #--- Save -----
    outbasedir = '/home/utsumi/temp/bams2020/ann-max'
    outdir = outbasedir + '/%s'%(vname)
    util.mk_dir(outdir)
    outpath= outdir + '/%s-%03d.npy'%(scen,ens)
    np.save(outpath, a3max)
    print outpath
