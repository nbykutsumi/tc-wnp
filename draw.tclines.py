# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import sys, os
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   detect_fsub import *
from   datetime import datetime, timedelta
import matplotlib.pyplot as plt
import config
import detect_func
import BestTrackTC
import util
import Cyclone
import ConstCyclone
import d4PDF
#--------------------------------------
#prj     = "HAPPI"
#model   = "MIROC5"
#run     = "C20-ALL-001"
#res     = "128x256"

prj     = "d4PDF"
model   = "__"
#scen    = 'HPB'
scen    = 'HPB'
ens     = 1
#run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
run     = "XX-%s-%03d"%(scen, ens)   # {expr}-{scen}-{ens}
res     = "320x640"
noleap  = False
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
wsDir     = wsbaseDir + '/%s'%(run)

iDTime = datetime(2004,1,1,0)
eDTime = datetime(2004,12,31,18)
lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(hours=6))

if prj=='d4PDF':
    a1lat   = d4PDF.Lat()
    a1lon   = d4PDF.Lon()
    ny      = len(a1lat)
    nx      = len(a1lon)


const  = ConstCyclone.Const(prj=prj, model=model)
thsst = 27
exrvort= 3.0*1e-5
tcrvort= 5.0*1e-5
thwcore= 3.0
thdura = 36
thwind = 15.
thwdif = -9999.

const['Lat']     = a1lat
const['Lon']     = a1lon
const['thsst']   = thsst + 273.15   # K
const['exrvort'] = exrvort
const['tcrvort'] = tcrvort
const['thwcore'] = thwcore
const['thdura']  = thdura
const['thwind']  = thwind
const['thwdif']  = thwdif

exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5
slabel = 'sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-wind-%02d-wdif-%d-du-%02d'%(thsst, exrvortout, tcrvortout, thwcore, thwind, thwdif, thdura)

cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#[[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,360]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[47,150]]

#----------------------------------

lonlatfontsize = 10.0
#lonrotation    = 90
lonrotation    = 0
miss_int= -9999

#------------------------
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
_, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
_, dtcppos= cy.mkInstDictC_objTC(iYM,eYM,varname='prepos')
#------------------------
print "Basemap"
figmap   = plt.figure(figsize=(7,7))
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

for DTime in lDTime:
  lxyz = dtcxy[DTime]
  lprepos = dtcppos[DTime]

  if len(lxyz)==0: continue

  for i,(x,y,z) in enumerate(lxyz):
    lon = a1lon[x]
    lat = a1lat[y]

    xpre,ypre = Cyclone.fortpos2pyxy(lprepos[i][2], nx, miss_int=-9999)
    if (xpre>0):
      lonpre,latpre = a1lon[xpre], a1lat[ypre]
    else:
      lonpre,latpre = lon, lat

    #if ((lllon<lon)&(lon<urlon)&(lllat<lat)&(lat<urlat)):
    #  print x,y, "***",lon, lat
    #  M.plot( lon, lat, "o")

    #------------------------------------
    lon1, lat1 = lon, lat
    lon2, lat2 = lonpre, latpre

    scol = 'r'
    if abs(lon1 - lon2) >= 180.0:
      #--------------
      if (lon1 > lon2):
        lon05_1  = 360.0
        lon05_2  = 0.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      elif (lon1 < lon2):
        lon05_1  = 0.0
        lon05_2  = 360.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      #--------------
      M.plot( (lon1, lon05_1), (lat1, lat05), linewidth=1, color=scol)
      M.plot( (lon05_2, lon2), (lat05, lat2), linewidth=1, color=scol)
      #--------------
    else:
      M.plot( (lon1, lon2), (lat1, lat2), linewidth=1, color=scol)



#-- meridians and parallels
gridflag = False
if gridflag == True:
  gridlinewidth = 1.0
elif gridflag == False:
  gridlinewidth = 0.0
meridians = 10.0
parallels = 10.0
M.drawmeridians(arange(0.0,360.0, meridians), labels=[0, 0, 0, 1], rotation=90, linewidth=gridlinewidth)
M.drawparallels(arange(-90.0,90.0, parallels), labels=[1, 0, 0, 0], linewidth=gridlinewidth)


#-- coastline ---------------
print "coastlines"
M.drawcoastlines(color="k")

#stitle  = '%s %03d %s-%s'%(scen, ens, iDTime,eDTime) + '\n' + slabel
stitle  = '%s %s-%s'%(run, iDTime,eDTime) + '\n' + slabel
plt.title(stitle)
figdir  = '/home/utsumi/temp/bams2020'
util.mk_dir(figdir)
figpath = figdir + '/temp.%s.png'%(run)
#figpath = figdir + '/temp-hk.png'
plt.savefig(figpath)
print figpath
plt.show()
sys.exit()



# %%

# %%
