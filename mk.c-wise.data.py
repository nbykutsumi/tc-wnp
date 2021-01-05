# %%
import matplotlib
#%matplotlib inline
#----------------------------------
import sys, os, pickle
#from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import myfunc.util as util
import calendar
from collections import deque
#import Cyclone
import socket
#--------------------------------------
#iY, eY = 1980,2010
iY, eY = 1990,2010
#iY, eY = 1990,1990
#iY, eY = 2001,2010
lYear = range(iY,eY+1)
lMon  = range(1,12+1)

#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lscen    = ["HPB","HPB_NAT"] # run={expr}-{scen}-{ens}
#lscen    = ["HPB"] # run={expr}-{scen}-{ens}
lscen    = ["HPB_NAT"] # run={expr}-{scen}-{ens}
#lens    = range(21,50+1)
lens    = range(3,50+1)
#lens    = range(1,2+1)
#lens    = range(38,50+1)
#lens    = range(41,50+1)
#lens    = range(1,50+1)
#lens    = [1]
#lens    = range(3,9+1)
res     = "320x640"
noleap  = False


detectName = 'wsd_d4pdf_20201209-py38'
#detectName = 'wsd_d4pdf_20200813'
#detectName = 'detect_20200401'
#config      = import_module("%s.config"%(detectName))
ConstCyclone= import_module("%s.ConstCyclone"%(detectName))
Cyclone     = import_module("%s.Cyclone"%(detectName))  # test
d4PDF       = import_module("%s.d4PDF"%(detectName))

hostname  = socket.gethostname()
if hostname=="shui":
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
elif hostname=="well":
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999

#************************************
# d4PDF (Objective detection)
#************************************
#lvar = ["age","dtlw","dtmd","dtup","dura","epos","idate","iedist","initland","initsst","ipos","land","lat","lon","nextpos","nowpos","prepos","slp","slp_mean_adj","slp_mean_box","sst","time","vortlw","vortlw_max_box","wmaxlw","wmaxup","x","y"]


lvar = ["dtlw","dtmd","dtup","lat","lon","slp","slp_mean_adj","slp_mean_box","vortlw","vortlw_max_box","wmaxlw","wmaxup","x","y"]

ddtype={
"dtlw": float32,
"dtmd": float32,
"dtup": float32,
"initland": float32,
"initsst": float32,
"lat": float32,
"lon": float32,
"slp": float32,
"slp_mean_adj": float32,
"slp_mean_box": float32,
"vortlw": float32,
"vortlw_max_box": float32,
"wmaxlw": float32,
"wmaxup": float32,
"dura":int32,
"x": int32,
"y": int32,
"tcmask":bool,
}
#lvar = lvar[:1]
#lMon = lMon[:1]
#lens = lens[:1]
#lscen= lscen[:1]

for scen in lscen:
    for ens in lens:

        try:
            del dout
        except NameError:
            pass

        dout = {var:{} for var in lvar}

        for Year in lYear:
            for Mon in lMon:
                print(ens,Year,Mon)
                clistdir = wsbaseDir + "/XX-%s-%03d/6hr/clist/%04d/%02d"%(scen,ens,Year,Mon)

                ipospath = clistdir + "/ipos.%04d.%02d.npy"%(Year,Mon)
                idatepath = clistdir + "/idate.%04d.%02d.npy"%(Year,Mon)
                iagepath  = clistdir + "/age.%04d.%02d.npy"%(Year,Mon)
                idurapath = clistdir + "/dura.%04d.%02d.npy"%(Year,Mon)
                iinitlandpath = clistdir + "/initland.%04d.%02d.npy"%(Year,Mon)
                iinitsstpath  = clistdir + "/initsst.%04d.%02d.npy"%(Year,Mon)

                aipos  = np.load(ipospath) 
                aidate = np.load(idatepath) 
                aage   = np.load(iagepath)
                adura  = np.load(idurapath)
                ailand = np.load(iinitlandpath)
                aisst  = np.load(iinitsstpath)


                dain = {}
                for var in lvar:  # test
                    ivarpath  = clistdir + "/%s.%04d.%02d.npy"%(var,Year,Mon)
                    dain[var] = np.load(ivarpath)

                nrec = len(dain[lvar[0]])

                for i in range(nrec):
                    ipos  = aipos[i] 
                    idate = aidate[i] 
                    age   = aage[i]
                    dura  = adura[i]

                    if dura < 36: continue
                    if dura > 720: continue  # longer than 30 days

                    try:
                        for var in lvar:
                            dat = dain[var][i]
                            dout[var][(idate,ipos)].append(dat)
                    except:
                        for var in lvar:
                            dat = dain[var][i]
                            dout[var][(idate,ipos)] = deque([dat])

                    if age==dura:
                        sidate = str(idate)
                        syear = sidate[:4]
                        smd   = sidate[4:8]
                        odir = wsbaseDir + "/XX-%s-%03d/6hr/cwise/%s/%s"%(scen,ens,syear,smd)
                        util.mk_dir(odir)

                        for var in lvar:
                            aout = np.asarray(dout[var][(idate,ipos)])
                            opath= odir + "/%s.%s.%s.bn"%(var,idate,ipos)
                            aout.astype(ddtype[var]).tofile(opath)
                            #np.save(opath, aout)
                        #-- single variables ---
                        iland = ailand[i]
                        isst  = aisst[i]

                        odurapath  = odir + "/dura.%s.%s.bn"%(idate,ipos)
                        oilandpath = odir + "/initland.%s.%s.bn"%(idate,ipos)
                        oisstpath  = odir + "/initsst.%s.%s.bn"%(idate,ipos)

                        dura.astype(ddtype["dura"]).tofile(odurapath)
                        iland.astype(ddtype["initland"]).tofile(oilandpath)
                        isst.astype(ddtype["initsst"]).tofile(oisstpath)
                        
                        #np.save(odurapath, dura)
                        #np.save(oilandpath, iland)
                        #np.save(oisstpath, isst)
                        #sys.exit()

                        del dout[var][(idate,ipos)]
                print(odir)


