from numpy import *
import d4PDF
import myfunc.util as util
import numpy as np

lYear = range(1980,2010+1)
lYear = [2010]
lscen = ["HPB"]
lens  = range(1,20+1)
lvname = ["TGEF"]

#dbbaseDir = "/home/utsumi/mnt/lab_work_hk03/d4PDF_GCM"
dbbaseDir = "/work/hk03/d4PDF_GCM"
d4mon = d4PDF.avr_mon_320x640(vtype="sfc", dbbaseDir=dbbaseDir)

for scen in lscen:
    for ens in lens:
        for vname in lvname:
            for Year in lYear:
                a2year = np.array([d4mon.load_ave_mon(vname=vname, scen=scen, ens=ens, Year=Year, Mon=Mon) for Mon in range(1,12+1)]).mean(axis=0)
                ny,nx = a2year.shape
                obasedir = "/home/utsumi/temp/bams2020/map-year"
                odir = obasedir + "/%s/%s.%03d"%(vname,scen,ens)
                util.mk_dir(odir)
                oname= odir + "/%s.%04d.%dx%d.npy"%(vname,Year,ny,nx)
                np.save(oname, a2year)
                print(oname)
        
