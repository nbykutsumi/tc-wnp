import numpy as np
from datetime import datetime, timedelta
import myfunc.util as util
import os, sys
import pickle

srcpath = '/home/utsumi/temp/bams2020/composite/obs/2019.10.NR.pickle'
with open(srcpath,'rb') as f:
    dat = pickle.load(f)

a1time = dat['time']
a1lat  = dat['lat']
a1lon  = dat['lon']
a1vort = dat['vort']
for i in range(len(a1time)):
    dtime= a1time[i]
    lat = a1lat[i]
    lon = a1lon[i]
    vort= a1vort[i]

    if ((dtime < datetime(2019,10,10,0))or(datetime(2019,10,13,18)<dtime)): continue
    if ((lat<0)or(lat>50)):continue
    if ((lon<130)or(150<lon)):continue
    print dtime,lon,lat,vort
