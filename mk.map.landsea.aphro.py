# %%
import numpy as np
import numpy.ma as ma
import APHRODITE
import matplotlib
import matplotlib.pyplot as plt
%matplotlib inline

#dbbaseDir = "/tank/utsumi/data/APHRO/APHRO_V1101"
dbbaseDir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101"

region= "MA"
res   = "050"
aph = APHRODITE.v1101(region=region, res=res, dbbaseDir=dbbaseDir)
ny  = aph.ny
nx  = aph.nx

Year = 2015
a3dat = aph.load_year(Year=Year)[0]
a2dat = ma.masked_less(a3dat,0).sum(axis=0)
a2landsea = np.ones([ny,nx],"int16")
a2landsea = ma.masked_where(a2dat.mask, a2landsea).filled(0)
#plt.imshow(a2landsea, origin="lower")
plt.imshow(ma.masked_less(a3dat[0],0), origin="lower")
plt.colorbar()
plt.show()

# %%
