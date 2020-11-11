import numpy as np
import d4PDF
import myfunc.util as util
import os, sys, pickle

lscen = ['HPB']
lens  = range(1,9+1)
compbasedir= '/home/utsumi/temp/bams2020/composite'
lvname  = ['prec']
a1latin = d4PDF.Lat()   # dlat ~ 0.5615674
a1lonin = d4PDF.Lon()  
lYear = range(2000,2010+1)
lMon  = range(1,12+1)

thsst = 27
thexrvort = 3 *1.0e-5
thtcrvort = 5 *1.0e-5
thwcore   = 3
thdura    = 36
thwind    = 15
thwdif    = -20

#--- Median of pmax ---
def load_a2th(vname):
    a3th = []
    a3th = np.concatenate([np.load('/home/utsumi/temp/bams2020/ann-max/%s/HPB-%03d.npy'%(vname,ens)) for ens in range(1,9+1)], axis=0)
    a2th = np.median(a3th, axis=0)
    return a2th

#----------------------
slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

for vname in lvname:
    a2th = load_a2th(vname)

    for scen in lscen:
        for ens in lens:
            for Year in lYear:
                for Mon in lMon:
                    compdir = compbasedir + '/%s/%s.%03d'%(slabel,scen,ens)
                    util.mk_dir(compdir)

                    comppath =  compdir + '/%04d.%02d.pickle'%(Year,Mon)

                    with open(comppath, 'rb') as f:
                        ddat = pickle.load(f)

                    a3var = ddat[vname]
                    a1lat = ddat['lat']
                    a1lon = ddat['lon']
                    a1dura= ddat['dura']
                    a1vort= ddat['vort']
                    a1tsum= ddat['tsum']
                    a1wup = ddat['wup']
                    a1wlw = ddat['wlw']

                    #-- Screening ----
                    a1flagdura = ma.masked_greater_equal(a1dura, thdura).mask

                    a1flagvort = ma.masked_greater_equal(a1vort, thtcrvort).mask

                    a1flagtsum = ma.masked_greater_equal(a1tsum, thwcore).mask

                    a1flagwind = ma.masked_greater_equal(a1wlw, thwind).mask

                    a1flagwdif = ma.masked_greater_equal(a1wlw - a1wup, thwdif).mask

                    a1flag = a1flagdura * a1flagvort * a1flagtsum * a1flagwind * a1flagwdif

                    if a1flag is np.bool_(False): continue

                    a1lat = a1lat[a1flag]
                    a1lon = a1lon[a1flag]
                    a2var = a3var[a1flag,:,:]

                    for i in range(len(a1lat)):
                        lat = a1lat[i]
                        lon = a1lon[i]


