# %%
import matplotlib
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
from numpy import *
from datetime import datetime, timedelta
import myfunc.util as util
import GSMaP


stpath = "/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.Ku.V06A.heightStormTop/2014/08/02/heightStormTop.1.002418.npy"

prpath = "/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.DPRGMI.V06A.9ave.surfPrecipTotRate/2014/08/02/surfPrecipTotRate.002418.npy"

st = ma.masked_less_equal(np.load(stpath), 0)
pr = ma.masked_less_equal(np.load(prpath), 0)[:,(103-83): (103-83)+(117-103+1)]

#plt.plot(st,pr,".")

for i in range(st.shape[0]):
    if ma.masked_less_equal(pr[i],0).count()>5:
        print(i)
# %%
