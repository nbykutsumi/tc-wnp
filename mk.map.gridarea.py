# %%
import matplotlib
#%matplotlib inline
import matplotlib.pyplot as plt
import d4PDF
from numpy import *
import math
import numpy as np

a1latbnd = d4PDF.LatBnd()
a1lonbnd = d4PDF.LonBnd()
ny = len(d4PDF.Lat())
nx = len(d4PDF.Lon())

a1lat0 = a1latbnd[:-1]
a1lat1 = a1latbnd[1:]
a1lon0 = a1lonbnd[:-1]
a1lon1 = a1lonbnd[1:]

Lon0, Lat0 = np.meshgrid(a1lon0, a1lat0)
Lon1, Lat1 = np.meshgrid(a1lon1, a1lat1)

pi = math.pi
R = 6378.136
e = 0.081819191
e2= 0.006694380

def calc_f(Lat):
    Lat = np.deg2rad(Lat)
    return 0.5*np.sin(Lat)/(1-e2*(np.sin(Lat)**2)) + 1/(4*e)*np.log(np.abs((1+e*np.sin(Lat))/(1-e*np.sin(Lat))))  # np.log: natural logarithm

f0 = calc_f(Lat0)
f1 = calc_f(Lat1)

a2area = pi*R**2*(1-e2)/180 * (f1 - f0) * np.abs(Lon1-Lon0)

opath = "./gridarea.%dx%d.npy"%(ny,nx)
np.save(opath, a2area.astype("float64"))
print(opath)


# %%
