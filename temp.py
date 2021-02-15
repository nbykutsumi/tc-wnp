# %%
import matplotlib
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import gzip
import struct
import calendar
import socket
import numpy.ma as ma
import APHRODITE

#dbbaseDir = "/tank/utsumi/data/APHRO/APHRO_V1101"
dbbaseDir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101"
Year = 1996
aphro = APHRODITE.v1101(region="MA", res="050", dbbaseDir=dbbaseDir)
aprec, aratio = aphro.load_year(Year)
print(aprec.shape)

#if calendar.isleap(Year):
#    nday = 366
#else:
#    nday = 365
#
#myhost = socket.gethostname()
#if myhost == "well":
#    srcdir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101/APHRO_MA/050deg"
#elif myhost=="shui":
#    srcdir = "/tank/utsumi/data/APHRO/APHRO_V1101/APHRO_MA/050deg"
#
#
#srcpath= srcdir + "/APHRO_MA_050deg_V1101.%04d.gz"%(Year)
#
#ny,nx = 140, 180
#with gzip.open(srcpath, "rb") as f:
#    #adat = np.fromfile(f, "float32").reshape(nday,2,ny,nx)
#    #adat = np.fromfile(f, "float32")
#    sdat = f.read()
#
#sfmt_all = "%df"%(ny*nx*2*nday)
#adat = np.array(struct.unpack(sfmt_all, sdat)).reshape(nday,2,ny,nx)
#print(adat[0,0])
plt.figure()
plt.imshow(ma.masked_less(aprec[100],0),origin="lower")
plt.colorbar()
plt.show
plt.figure()
plt.imshow(ma.masked_less(aratio[100],0),origin="lower")
plt.colorbar()
plt.show

#plt.figure()
#plt.imshow(ma.masked_less(adat[50,0],0),origin="lower")
#plt.colorbar()
#plt.show
#plt.figure()
#plt.imshow(ma.masked_less(adat[50,1],0),origin="lower")
#plt.colorbar()
#plt.show






# %%
