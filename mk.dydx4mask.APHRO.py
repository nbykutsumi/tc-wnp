# %%
import numpy as np
from math import sin, cos, acos
import myfunc.util as util
from myfunc.util import  calc_dist_gc
import os, sys
import APHRODITE

region="MA"
res   ="050"
aph = APHRODITE.v1101(region=region, res=res)
ny   = aph.ny
nx   = aph.nx
a1lat= aph.Lat

dlat = {"050":0.5, "025":0.25}[res]
dlon = {"050":0.5, "025":0.25}[res]

thradkm= 500.
#thradkm= 200.
dykm = calc_dist_gc(0, dlat, 0, 0)
ndy = int(thradkm/dykm)

otable = np.full([ny,2*ndy+1],-9999, "int32")
print(ny,ndy)

for i, lat0 in enumerate(a1lat):
    print(lat0)
    for j,djy in enumerate(range(-ndy,ndy+1)):
        if ((i+djy)<0)or((i+djy)>(ny-1)): continue

        lat1 = a1lat[i+djy]

        for dkx in range(int(nx/2),0-1,-1):
            if (djy==0)&(dkx==0): continue
            lon1=dlon*dkx
            #print(lat0,lat1,0,lon1,dkx)
            dist = calc_dist_gc(lat0,lat1,0,lon1)
            #print(lat0,j,lat1,dkx,lon1,dist)
            if dist<thradkm:
                otable[i,j]=dkx
                break

outpath = "./tab.dydx4mask.APHRO.%s.%sdeg.nrady-%03d.%04dkm.npy"%(region,res,ndy,thradkm)
np.save(outpath, otable.astype("int32"))
print(outpath)
print(otable[0])
print(otable[100])
# %%
