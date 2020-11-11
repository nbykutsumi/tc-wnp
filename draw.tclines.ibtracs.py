# %%
import matplotlib
%matplotlib inline
from numpy import *
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import IBTrACS
import util
from datetime import datetime, timedelta
import os, sys
#-------------------------------
gridflag = False
idtime = datetime(2004,1,1,0)
edtime = datetime(2004,12,31,18)
#idtime = datetime(2000,1,1,0)
#edtime = datetime(2000,12,31,18)

#[[lllat,lllon],[urlat,urlon]] = [[30,120],[47,150]]
[[lllat,lllon],[urlat,urlon]] = [[0,100],[47,150]]
ib = IBTrACS.IBTrACS()
dtc = ib.ret_dlonlat(idtime,edtime, dhours=6, lvaridx=[0])
ldtime= dtc.keys()
ltc = dtc.values()

ltc = zip(ldtime, ltc)
ltc = sorted(ltc, key=lambda x:x[0]) # sort by datetime
ltc = zip(*ltc)[1]  # take only values
ltc = [x for l in ltc for x in l]

a1lon, a1lat, a1tcid = map(np.array, zip(*ltc))
#print a1tcid
#---- Screening ----------
a1flaglat = ma.masked_inside(a1lat, lllat-5, urlat+5).mask
a1flaglon = ma.masked_inside(a1lon, lllon-5, urlon+5).mask
a1flag    = a1flaglat * a1flaglon

a1lon = a1lon[a1flag]
a1lat = a1lat[a1flag]
a1tcid= a1tcid[a1flag]

ddat = {}
for (lon,lat,tcid) in zip(a1lon, a1lat, a1tcid):
    if tcid not in ddat.keys():
        ddat[tcid] = [[lon,lat]]
    else:
        ddat[tcid].append([lon,lat])
  
  
#**********************************************
#---- Draw tc track lines -------------------
#**********************************************
#------------------------
# Basemap
#------------------------
figmap   = plt.figure(figsize=(7,7))
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.7])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
#---- para ---------------
#scol = "r"
scol = "0.1"
#*************\***********
for tcid in ddat.keys():
    a1lon, a1lat = zip(*ddat[tcid])

    M.plot(a1lon, a1lat, '-')


#-- coastline ---------------
print "coastlines"
M.drawcoastlines(color="k")

#-- meridians and parallels
if gridflag == True:
  gridlinewidth = 1.0
elif gridflag == False:
  gridlinewidth = 0.0
meridians = 10.0
parallels = 10.0
M.drawmeridians(arange(0.0,360.0, meridians), labels=[0, 0, 0, 1], rotation=90, linewidth=gridlinewidth)
M.drawparallels(arange(-90.0,90.0, parallels), labels=[1, 0, 0, 0], linewidth=gridlinewidth)

stitle = 'Best track' + ' %s-%s'%(idtime,edtime) + '\n'+ '[[%d,%d],[%d,%d]]'%(lllat,lllon,urlat,urlon)
plt.title(stitle)
plt.show()

    ##----------
    #if i == len(lines)-1:
    #  continue
    ##----------
    #line_now = lines[i]
    #line_nxt = lines[i+1]
    #lat1     = line_now[0]
    #lon1     = line_now[1] 
    #lat2     = line_nxt[0]
    #lon2     = line_nxt[1] 

    ##------------------------------------
    #if abs(lon1 - lon2) >= 180.0:
    #  #--------------
    #  print tcid,lat1,lat2, lon1,lon2


    #  if (lon1 > lon2):
    #    lon05_1  = 360.0
    #    lon05_2  = 0.0
    #    lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
    #  elif (lon1 < lon2):
    #    lon05_1  = 0.0
    #    lon05_2  = 360.0
    #    lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)

    #  #--------------
    #  M.plot( (lon1, lon05_1), (lat1, lat05), linewidth=1, color=scol)
    #  M.plot( (lon05_2, lon2), (lat05, lat2), linewidth=1, color=scol)
    #  #--------------
    #else:
    #  M.plot( (lon1, lon2), (lat1, lat2), linewidth=1, color=scol)

##-- coastline ---------------
#print "coastlines"
#M.drawcoastlines(color="k")
#
##-- meridians and parallels
#if gridflag == True:
#  gridlinewidth = 1.0
#elif gridflag == False:
#  gridlinewidth = 0.0
#meridians = 10.0
#parallels = 10.0
#M.drawmeridians(arange(0.0,360.0, meridians), labels=[0, 0, 0, 1], rotation=90, linewidth=gridlinewidth)
#M.drawparallels(arange(-90.0,90.0, parallels), labels=[1, 0, 0, 0], linewidth=gridlinewidth)
#
##---- title --
#stitle  = "best track\n" + "%04d-%04d season:%s"%(iyear, eyear, season)
#axmap.set_title(stitle)
##--- save ---
#odir   = odir_root + "/%04d-%04d.%s"%(iyear,eyear, season)
#ctrack_func.mk_dir(odir) 
#soname = odir + "/tclines.ibtracs_all.v03r04.%04d-%04d.%s.png"%(iyear,eyear,season)
#soname = odir + "/tclines.ibtracs_all.v03r04.%04d-%04d.%s.png"%(iyear,eyear,season)
#plt.savefig(soname)
#print soname
##--------------------



# %%
