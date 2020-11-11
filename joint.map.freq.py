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
iY, eY = 2000,2010

#lthsst  = [27,28]
lthsst  = [27,27.5,28]
lexrvort = np.array([3])*1.0e-5
#ltcrvort = np.array([5])*1.0e-5
ltcrvort = np.array([3])*1.0e-5
lthwcore= [0]
lthdura = [36]
lthwind = [10,11,12, 13]
#lthwdif = [-9999]
#lthwdif = [-9999]
lthwdif = [-9999, -30,-20,-10]



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

i = -1
ddat = {}
for key in lkey:
    i=i+1

    [thsst,exrvort,tcrvort,thwcore,thdura,thwind,thwdif] = key

    #slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)
    slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst*10, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

    figdir = '/home/utsumi/temp/bams2020/fig/map-freq'
    figpath = figdir + '/map.freq.tc.obj.%s.%04d-%04d.ave.png'%(slabel, iY, eY)


    try:
        iimg = Image.open(figpath)
        a2array = asarray(iimg)[iy:ey, ix:ex]
    except IOError:
        a2array = a2array*0 + 255
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
outpath = figdir + '/joint.map.freq.tc.obj.%04d-%04d.ave.png'%(iY, eY)
oimg.save(outpath)
print outpath




# %%
