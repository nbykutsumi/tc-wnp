# %%
#import matplotlib
#import matplotlib.pyplot as plt
#matplotlib.use("Agg")
#%matplotlib inline
import numpy as np
import scipy
import d4PDF
import sys, os
import myfunc.util as util
from datetime import datetime, timedelta
import socket
from importlib import import_module
from numpy import ma
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
#lscen   = ['HPB']  # run={expr}-{scen}-{ens}
lscen   = ['HPB_NAT'] # run={expr}-{scen}-{ens}
#lscen   = ['HPB','HPB_NAT'] # run={expr}-{scen}-{ens}

lens    = list(range(1,50+1))
#lens    = list(range(21,50+1))
#lens    = [1]
miss_out = -9999.

iY,eY = 1990,2010
#iY,eY = 1990,1990
lY  = range(iY,eY+1)
ny,nx = 320, 640
hostname = socket.gethostname()
if hostname =='shui':
    wsbaseDir = '/tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/work/hk03/d4PDF_GCM'
elif hostname=='well':
    wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'
    d4pdfdir = '/home/utsumi/mnt/lab_work/hk03/d4PDF_GCM'

detectName = 'wsd_d4pdf_20201209-py38'
d4PDF       = import_module("%s.d4PDF"%(detectName))
d4sfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=d4pdfdir)

a2land = d4PDF.load_topo_TL319(dbbaseDir=d4pdfdir, vname='ratiol', miss_fill=-9999.)

for scen in lscen:
    for ens in lens:
        for Year in lY:

            a3prec_day = ma.concatenate([d4sfc.load_6hr_mon("PRECIPI", scen, ens, Year, Mon).reshape(-1,4,ny,nx).mean(axis=1) for Mon in range(1,12+1)], axis=0)

            a2max = a3prec_day.max(axis=0).filled(miss_out)   # No unit change. kg/m2/sec
            #--- Save -----
            maxdir = "/home/utsumi/temp/bams2020/max-prec-d/%s.%03d"%(scen,ens)
            maxpath= maxdir + "/ann-max-prec.%04d.npy"%(Year)
            util.mk_dir(maxdir)
            np.save(maxpath, a2max)
            print(maxpath)
# %%
