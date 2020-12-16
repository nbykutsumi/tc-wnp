# %%
import matplotlib
import numpy as np
from numpy import ma
from datetime import datetime, timedelta
import calendar
from d4PDF import avr_mon_320x640
import matplotlib.pyplot as plt
#%matplotlib inline

myhost = 'shui'
#myhost = 'well'
if myhost == 'well':
    hisbasedir = '/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM'
    natbasedir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
elif myhost == 'shui':
    hisbasedir = '/work/hk03/d4PDF_GCM'
    natbasedir = '/work/hk01/d4PDF_GCM'

lyear = range(1980,2010+1)
#lyear = range(1999,2010=1)

vname = 'TA'
lens = range(1,10+1)

d4 = avr_mon_320x640(vtype='sfc', dbbaseDir=hisbasedir)

lregion = ['WNP', 'NEA']
region = 'WNP'
dbbox = {
    'WNP':[[20,120],[47,150]],
    'NEA':[[30,120],[47,150]],
    }

[[lllat,lllon],[urlat,urlon]] = dbbox[region]

a1lat = d4.Lat
a1lon = d4.Lon
a2lon, a2lat = np.meshgrid(a1lon, a1lat)
a2lonmask = ma.masked_outside(a2lon, lllon, urlon).mask
a2latmask = ma.masked_outside(a2lat, lllat, urlat).mask

a2regionmask = a2lonmask + a2latmask

lhis = []
scen = 'HPB'
for ens in lens:
    ltmp = []
    for year in lyear:
        print scen,ens, year
        a2ta = np.array([d4.load_ave_mon(vname, scen, ens, year,mon) for mon in range(1,12+1)]).mean(axis=0)
        ta   = ma.masked_where(a2regionmask, a2ta).mean() - 273.15
        #ta = 18 + 2*np.sin((0.5*year*(ens*0.2)*np.pi))  # test
        ltmp.append(ta)
    lhis.append(ltmp)
a2his = np.array(lhis)


d4 = avr_mon_320x640(vtype='sfc', dbbaseDir=natbasedir)
lnat = []
scen = 'HPB_NAT'
for ens in lens:
    ltmp = []
    #for year in lyear:
    for year in lyear:
        print scen, ens, year
        if year >=2000:
            a2ta = np.array([d4.load_ave_mon(vname, scen, ens, year,mon) for mon in range(1,12+1)]).mean(axis=0)
            ta   = ma.masked_where(a2regionmask, a2ta).mean() - 273.15
            #ta = 18 + 2*np.sin((year*(ens*0.2)*np.pi+1))  # test
        else:
            ta = np.nan
        ltmp.append(ta)
    lnat.append(ltmp)
a2nat = np.array(lnat)
#--------------------------------
# Figure
#--------------------------------
fig = plt.figure(figsize=(8,5))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
a1year = np.array(lyear)

for iens in range(len(lens)):
    ax.plot(a1year, a2his[iens], color='0.5', linestyle='-', linewidth=0.5)
    ax.plot(a1year, a2nat[iens], color='b', linestyle='-', linewidth=0.5)

ax.plot(a1year, a2his.mean(axis=0), color='k', linestyle='-', linewidth=2)
ax.plot(a1year, a2nat.mean(axis=0), color='b', linestyle='-', linewidth=2)

stitle = '2m temperature [deg. C] %s'%(region)
plt.title(stitle)
print a2nat.shape
print a2his.shape
print len(lyear)
print a2his
sdomain = 'lat.%03d-%03d.lon.%04d-%04d'%(lllat,urlat,lllon,urlon)
figdir = '/home/utsumi/temp/bams2020/fig/%s'%(sdomain)
figpath = figdir + '/ta.ts.png'
plt.savefig(figpath)
print figpath
plt.show()
# %%
