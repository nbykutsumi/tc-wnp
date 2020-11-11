# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
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
lens    = range(1,9+1)
#lens    = range(2,9+1)
#lens    = [1]

#iYM = [2000,1]
#eYM = [2000,12]

iYM = [2001,1]
eYM = [2010,12]
lYM = util.ret_lYM(iYM, eYM)

lscen   = ['HPB','HPB_NAT']
#lscen   = ['HPB_NAT']
#lscen   = ['HPB']

d4pdfdir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)
d4atm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir=d4pdfdir)
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

y0 = bisect_left(a1latin, lllat)
y1 = bisect_left(a1latin, urlat)
x0 = bisect_left(a1lonin, lllon)
x1 = bisect_left(a1lonin, urlon)

dbins = {
    'prec': np.arange(0,100+0.01,2.5),
    'wind': np.arange(0,100+0.1, 2),
    }
print 'y0,y1',y0,y1
print 'x0,x1',x0,x1

for scen in lscen:
    for ens in lens:
        for (Year,Mon) in lYM:
            #-- load d4PDF variables -----
            a2prec = d4sfc.load_6hr_mon(vname='PRECIPI', scen=scen, ens=ens, Year=Year, Mon=Mon)[:,y0:y1+1,x0:x1+1] *60*60
            a2u    = d4sfc.load_6hr_mon(vname='UAOPN', scen=scen, ens=ens, Year=Year, Mon=Mon) [:,y0:y1+1,x0:x1+1]
            a2v    = d4sfc.load_6hr_mon(vname='VAOPN', scen=scen, ens=ens, Year=Year, Mon=Mon) [:,y0:y1+1,x0:x1+1]
            a2wind = np.sqrt(a2u**2 + a2v**2)

            #-- Histogram -----
            a2hist_prec, _ = np.histogram(a2prec, bins=dbins['prec'])
            a2hist_wind, _ = np.histogram(a2wind, bins=dbins['wind'])

            ddat = {}
            ddat['prec'] = a2hist_prec
            ddat['wind'] = a2hist_wind
            ddat['bins-prec'] = dbins['prec']
            ddat['bins-wind'] = dbins['wind']
            ddat['BBox'] = [[lllat,lllon],[urlat,urlon]]

            outdir = '/home/utsumi/temp/bams2020/histogram/lat.%03d-%03d.lon.%04d-%04d/gen.%s.%03d'%(lllat,urlat,lllon,urlon, scen, ens)
            util.mk_dir(outdir)

            outpath = outdir + '/%04d.%02d.pickle'%(Year,Mon)
            with open(outpath, 'wb') as f:
                pickle.dump(ddat, f)
            print outpath



# %%
