# %%
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline

prec = np.load('/home/utsumi/mnt/lab_tank/utsumi/hometemp/bams2020/composite/obs-era5-tc/1980/prec.1980.07.npy')
wind = np.load('/home/utsumi/mnt/lab_tank/utsumi/hometemp/bams2020/composite/obs-era5-tc/1980/wind.1980.07.npy')

a2prec=prec.mean(axis=0)
a2wind=wind.mean(axis=0)
plt.imshow(a2prec,origin='lower')
plt.plot(9,9,'x')
plt.colorbar()
plt.show()

plt.imshow(a2wind,origin='lower')
plt.plot(9,9,'x')
plt.colorbar()
plt.show()

# %%
