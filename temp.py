# %%
import matplotlib.pyplot as plt
#%matplotlib inline

import numpy as np
from math import sin, cos, acos
import myfunc.util as util
from myfunc.util import  calc_dist_gc
import os, sys
from detect_fsub import *

ny,nx = 320, 640
a1xpy = np.array([200, 200,250])
a1ypy = np.array([20,  100,150])
a2table = np.load("./tab.dydx4mask.d4PDF.nrady-008.0500km.npy")

a2out = detect_fsub.mk_a2mask_with_table(a2table.T, a1xpy, a1ypy, nx, ny).T

plt.imshow(a2out, origin="lower")
plt.show()
print(a2out.shape)
print(a2out)
print(a2out.max())

# %%
