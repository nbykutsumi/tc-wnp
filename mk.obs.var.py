# %%
#%matplotlib inline
import numpy as np
from numpy import ma
import myfunc.util as util
import myfunc.grids as grids
import Cyclone
import ConstCyclone
import IBTrACS, JRA55
from datetime import datetime, timedelta
import sys, os, calendar
import pickle
#************************************************
#obsflag = True
obsflag = False
#simflag = True
simflag = False


#************************************************
# Functions
#------------------------------------------------
#************************************************
# IBTrACS
#------------------------------------------------
#iYM  = [1980,1]
iYM  = [2007,4]
eYM  = [2019,12]
lYM  = util.ret_lYM(iYM,eYM)

bst = IBTrACS.IBTrACS()


jrabaseDir = '/home/utsumi/mnt/lab_data2/JRA55/'
jra = JRA55.anl_p125(dbbaseDir=jrabaseDir)
vname = 'relv'
a1latjra = JRA55.Lat125(crd='sa')
a1lonjra = JRA55.Lon125(crd='sa')
dlat,dlon= JRA55.dlatlon125()
latjra0  = a1latjra[0]
lonjra0  = a1lonjra[0]

for (Year,Mon) in lYM:
    eDay = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1,0)
    eDTime = datetime(Year,Mon,eDay,18)
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=6))
    bstlonlat = bst.ret_dlonlat(iDTime, eDTime)

    dout = {}
    for DTime in lDTime:
        print DTime
        nyjra = len(JRA55.Lat125())
        avar = jra.load_6hr(vname=vname, DTime=DTime, plev=850, crd='sa')
        avar[:nyjra/2] = (-ma.masked_equal(avar[:nyjra/2], -9999.)).filled(-9999.) # Flip signs of southern hemisphere
        avar = grids.karnel_pooling_map2D_global(avar, dy=1,dx=1,func='max',miss_in=-9999.)  # Find maximum in 3x3 grid boxes
        llonlat = bstlonlat[DTime]
        lout = []
        for (lon,lat) in llonlat:
            if lon<0: lon=lon+360
            y = int((lat - latjra0)/dlat)
            x = int((lon - lonjra0)/dlon)

            lout.append(avar[y,x])

        dout[DTime] = lout

    #-- Save -------
    pickledir = '/home/utsumi/temp/bams2020/tc-obs-var/%s'%(vname)
    util.mk_dir(pickledir)
    picklepath = pickledir + '/%04d.%02d.pickle'%(Year,Mon)
    with open(picklepath, 'wb') as f:
        pickle.dump(dout, f)
    print picklepath

