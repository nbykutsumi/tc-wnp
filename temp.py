# %%
import matplotlib.pyplot as plt
%matplotlib inline
import numpy as np
from numpy import *
import d4PDF
from bisect import bisect_left
import sys

srcpath = "/home/utsumi/temp/bams2020/tc-prec-0500km/HPB.001/timesper-1day.1990.npy"
a3num=np.load(srcpath)
for i in range(a3num.shape[0]):
    a2num = a3num[i]
    if a2num.max()==4:
        plt.imshow(a2num,origin="lower")
        plt.colorbar()
        plt.show()
        sys.exit()
# %%
