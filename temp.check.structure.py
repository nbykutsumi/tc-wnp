# %%
import matplotlib
matplotlib.use('Agg')
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
#--------------------------------------
#-----------------
iYM     = [2000,1]
eYM     = [2009,12]
lYM  = util.ret_lYM(iYM, eYM)
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
scen    = 'HPB' # run={expr}-{scen}-{ens}
#lens    = range(1,9+1)
ens    =  1
res     = "320x640"
noleap  = False


detectName = 'wsd_d4pdf_20200428'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))
d4PDF       = import_module("%s.d4PDF"%(detectName))
IBTrACS     = import_module("%s.IBTrACS"%(detectName))

wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
outdir  = '/home/utsumi/temp/bams2020'
#----------------------------------
miss_int= -9999
[[lllat,lllon],[urlat,urlon]] = [[0,100],[60,180]]

a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)
#************************************

thsst = 26 # degC
thsstK = thsst + 273.15
exrvort   = 3*1.0e-5
tcrvort   = 3*1.0e-5
thwcore   = -1
thdura    = 0

#************************************
def solve_time(stime):
  year = int( stime/10**6 )
  mon  = int( (stime - year*10**6)/10**4 )
  day  = int( (stime - year*10**6 - mon*10**4)/10**2)
  hour = int( (stime - year*10**6 - mon*10**4 - day*10**2) )
  return year, mon, day, hour

def fortpos2pyxy(number, nx, miss_int):
  if (number == miss_int):
    iy_py = miss_int
    ix_py = miss_int
  else:
    iy_py = int((number-1.0)/nx)      # iy_py = 0,1,2,..
    ix_py = int(number - nx*iy_py-1)   # ix_py = 0,1,2,..
  #----
  return ix_py, iy_py
#***************************************************

wsDir= wsbaseDir + '/%s-%s-%03d'%(expr, scen, ens)
const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = d4PDF.Lat()
const['Lon'] = d4PDF.Lon()
cy  = Cyclone.Cyclone(baseDir=wsDir, const=const)

exrvortout = exrvort*1.0e+5
tcrvortout = tcrvort*1.0e+5

dictExC   = {}
dictTC    = {}
da1       = {}

lvar      = ["dura","nowpos","time","vortlw","dtlow","dtmid","dtup","initsst","initland","wmeanup","wmeanlow"]  # wmean added.

lout = []
for Year,Mon in lYM:
    print Year,Mon
    for var in lvar:
        da1[var]  = cy.load_clist(var, Year, Mon)


    nlist    = len(da1["dura"])
    for i in range(nlist):
        dura        = da1["dura"    ][i]
        nowpos      = da1["nowpos"  ][i]
        time        = da1["time"    ][i]
        #iedist      = da1["iedist"  ][i]
        #rvort       = abs(da1["vortlw"   ][i])
        rvort       = da1["vortlw"   ][i]
        dtlow       = da1["dtlow"   ][i]
        dtmid       = da1["dtmid"   ][i]
        dtup        = da1["dtup"    ][i]
        initsst     = da1["initsst" ][i]
        initland    = da1["initland"][i]
        #nextpos     = da1["nextpos" ][i]
        wup         = da1["wmeanup" ][i]
        wlow        = da1["wmeanlow"][i]

        #---- dura -------
        if dura < thdura:
          #print "dura",dura,"<",thdura
          continue

        #---- time -------
        Year,Mon,Day,Hour = solve_time(time)
        DTime             = datetime(Year,Mon,Day,Hour)

        #---- nowpos  ----
        x,y               = fortpos2pyxy(nowpos, nxin, -9999)

        lon = a1lonin[x]
        lat = a1latin[y]

        if (lon<lllon)or(urlon<lon):
            continue
        if (lat<lllat)or(urlat<lat):
            continue

        #---- thrvort ----
        if rvort < tcrvort:
          continue

        #---- thwcore ----
        if dtup + dtmid + dtlow < thwcore:
        #if (dtup <thwcore)or(dtmid<thwcore)or(dtmid<thwcore):
          #print "thwcore",dtup+dtmid+dtlow,"<",thwcore
          continue
        dtsum = dtup+dtmid+dtlow
        #---- wup & wlow --
        #if wup > wlow:
        #  try:
        #    dictExC[DTime].append(oList)
        #  except KeyError:
        #    dictExC[DTime] = [oList]

        #  #print "wup > wlow !!"
        #  continue

        #---- initsst ----
        if initsst < thsstK:
          #print "initsst",initsst,"<",thinitsst
          continue

        initsstC= initsst - 273.15
        #---- initland ----
        if initland > 0.0:
          #print "initland",initland,">",0.0
          continue

        #--------------------
        ltmp = [lat,lon,Mon,Day,initsstC,dura,dtsum,dtlow,dtmid,dtup,rvort,wlow,wup]

        lout.append(ltmp)

label= ['lat','lon','mon','day','initsst','dura','dtsum','dtlow','dtmid','dtup','rvort','wlow','wup']

lout.insert(0, label)
sout = util.list2csv(lout)

slabel = '%s-%s-sst-%d.ex-%.2f.tc-%.2f.wc-%.1f-du-%02d-en-%03d'%(expr, scen, thsst, exrvortout, tcrvortout, thwcore, thdura, ens)

csvpath = outdir + '/table.%s.csv'%(slabel)
f=open(csvpath,'w'); f.write(sout); f.close()
print csvpath

