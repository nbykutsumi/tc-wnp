# %%
import matplotlib
#%matplotlib inline

from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket
import numpy as np
#*******************************

iy =1 # top
ey = -1
ix =10
ex = -30
#**********
target  = 'point'
lregion = ['WNP', 'NEA']
dbbox = {
    'WNP':[[20,120],[47,150]],
    'NEA':[[30,120],[47,150]],
    }

#lthsst  = [27,28]
lthsst  = [27,27.5,28]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
lthdura = [36]
lthwind = [10,11,12,13]
#lthwdif = [-9999,0]
lthwdif = [-9999,-30,-20,-10]

for region in lregion:
    [[lllat,lllon],[urlat,urlon]] = dbbox[region]

    lkey = [[thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif]
            for thsst   in lthsst
            for exrvort in lexrvort
            for tcrvort in ltcrvort
            for thwcore in lthwcore
            for thdura in lthdura
            for thwind in lthwind
            for thwdif in lthwdif
            ]
    exrvortout = exrvort*1.0e+5
    tcrvortout = tcrvort*1.0e+5

    print len(lkey)
    #------------------------------
    # Monthly
    #------------------------------
    i = -1
    ddat = {}
    for key in lkey:
        i=i+1

        [thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif] = key
    
        #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

        sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)

        figdir = '/home/utsumi/temp/bams2020/fig/%s'%(sdomain)
        figpath = figdir + '/count.mon.%s.%s.%s.png'%(target, slabel, region)

        iimg = Image.open(figpath)

        a2array = asarray(iimg)[iy:ey, ix:ex]
        ddat[i] = a2array

    ddat[-1] = a2array *0 + 255

    nx = 8
    ny = 6
    a2oarray = []
    i = -1
    for y in range(ny):
        a2line = []
        for x in range(nx):
            i = i+1
            if i<len(lkey):
                a2line.append(ddat[i])
            else:
                a2line.append(ddat[-1])
        a2line = np.concatenate(a2line, axis=1)
        a2oarray.append(a2line)
    a2oarray = np.concatenate(a2oarray, axis=0)
    #a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3],ddat[4],ddat[5]],axis=1)
    #a2line1 = concatenate([ddat[6],ddat[7],ddat[8],ddat[9],ddat[10],ddat[11]],axis=1)
    #a2line2 = concatenate([ddat[12],ddat[13],ddat[14],ddat[15],ddat[16],ddat[17]],axis=1)
    #a2line3 = concatenate([ddat[18],ddat[19],ddat[20],ddat[21],ddat[22],ddat[23]],axis=1)
    #a2line4 = concatenate([ddat[24],ddat[25],ddat[26],ddat[27],ddat[28],ddat[29]],axis=1)
    #a2line5 = concatenate([ddat[30],ddat[31],ddat[32],ddat[33],ddat[34],ddat[35]],axis=1)
    #a2oarray = np.concatenate([a2line0, a2line1, a2line2, a2line3,a2line4,a2line5], axis=0)

    oimg    = Image.fromarray(a2oarray)
    outpath = figdir + '/joint.count.mon.%s.%s.png'%(target, region)
    oimg.save(outpath)
    print outpath

    #------------------------------
    # Annual
    #------------------------------
    i = -1
    ddat = {}
    for key in lkey:
        i=i+1

        [thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif] = key
    
        #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
        slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

        sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)

        figdir = '/home/utsumi/temp/bams2020/fig/%s'%(sdomain)
        figpath = figdir + '/count.year.%s.%s.%s.png'%(target, slabel, region)

        iimg = Image.open(figpath)

        a2array = asarray(iimg)[iy:ey, ix:ex]
        ddat[i] = a2array

    ddat[-1] = a2array *0 + 255

    nx = 8
    ny = 6
    a2oarray = []
    i = -1
    for y in range(ny):
        a2line = []
        for x in range(nx):
            i = i+1
            if i<len(lkey):
                a2line.append(ddat[i])
            else:
                a2line.append(ddat[-1])
        a2line = np.concatenate(a2line, axis=1)
        a2oarray.append(a2line)
    a2oarray = np.concatenate(a2oarray, axis=0)

    oimg    = Image.fromarray(a2oarray)
    outpath = figdir + '/joint.count.year.%s.%s.png'%(target, region)
    oimg.save(outpath)
    print outpath

    
    

# %%
   
    

# %%
