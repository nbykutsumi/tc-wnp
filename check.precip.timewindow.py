import numpy as np
from datetime import datetime, timedelta
import pygrib
import d4PDF
import util

daypath = "/home/utsumi/temp/bams2020/temp/sfc_avr_day_HPB_m001_199901.grib"
year = 1999
mon  = 1
day=1
#y,x = 262, 14
#y,x = 276, 44
y,x = 264, 349

with pygrib.open(daypath) as grbs:
    grb = grbs.message(day*2)
    #a2day = grbs.select(indicatorOfParameter="RI")[0].values
    a2day = grb.values*60*60*24
#print(a2day)
#print(a2day.shape)
#a2daymm = a2day*60*60*24
#for y in range(320):
#    for x in range(640):
#        if a2daymm[y,x] >5:
#            print(y,x, a2daymm[y,x])



six    = d4PDF.snp_6hr_2byte(vtype="sfc")

ldtimef = util.ret_lDTime(datetime(year,mon,day,0),datetime(year,mon,day,18),timedelta(hours=6))
a3sixf = np.array([six.load_6hr(vname="PRECIPI", scen="HPB", ens=1, DTime=DTime) for DTime in ldtimef])
a2sixf = a3sixf.sum(axis=0)*60*60*6

ldtimeb = util.ret_lDTime(datetime(year,mon,day,6),datetime(year,mon,day+1,0),timedelta(hours=6))
a3sixb = np.array([six.load_6hr(vname="PRECIPI", scen="HPB", ens=1, DTime=DTime) for DTime in ldtimeb])
a2sixb = a3sixb.sum(axis=0)*60*60*6


print(a2day[y,x])
print(a2sixf[y,x])
print(a2sixb[y,x])
print(a2day)
print(a2sixf)
print(a2sixb)

