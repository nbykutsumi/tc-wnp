import shutil
import os
import glob

lscen = ["HPB","HPB_NAT"]
#lscen = ["HPB_NAT"]
#lens = range(1,50+1)
lens = range(3,50+1)
for ens in lens:
    for scen in lscen:
        ssearch = "/tank/utsumi/WS/d4PDF_GCM/XX-%s-%03d*/6hr/cwise"%(scen,ens)
        lsrcpath = sorted(glob.glob(ssearch))
        
        for srcpath in lsrcpath:
            print(srcpath)
            shutil.rmtree(srcpath)

