# %%
import numpy as np
from numpy import ma
import regionmask
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib
import APHRODITE
import socket
import sys, os
matplotlib.use("Agg")
%matplotlib inline

hostname=socket.gethostname()
if hostname=="shui":
    dbbaseDir = "/tank/utsumi/data/APHRO/APHRO_V1101"
elif hostname=="well":
    dbbaseDir = "/home/utsumi/mnt/lab_tank/utsumi/data/APHRO/APHRO_V1101"
else:
    print("check hostname",hostname)
    sys.exit()

aphro = APHRODITE.v1101(region="MA", res="050", dbbaseDir=dbbaseDir)

a1latcnt = aphro.Lat
a1loncnt = aphro.Lon
a1latbnd = aphro.LatBnd
a1lonbnd = aphro.LonBnd

[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,150]]   # for figure

#JP = np.array([[123,25],[123,23],[130,23],[135,30],[142,30],[150,46],[140,46],[135,40],[129,34]])   # [lon, lat]
#SJ = np.array([[123,25],[123,23],[130,23],[135,30],[142,30],[150,46],[140,46],[135,40],[129,34]])   # [lon, lat]
NJ = np.array([[135,38],[150,38],[150,46],[138,46]])
SJ = np.array([[135,38],[123,30],[123,21],[133,21],[133,30],[143,30],[143,38]])
KR = np.array([[123,40],[123,30],[130,35],[130,40]])
EC = np.array([[116,42],[116,26],[123,26],[123,42]])
SC = np.array([[107,23],[107,18],[112,18],[123,21],[123,26],[118,26]])
IC = np.array([[105,23],[105,18],[107,12],[105,10],[105,7],[112,12],[112,18],[107,18],[107,23]])
PH = np.array([[119,20],[119,15],[124,5],[128,5],[128,10],[123,20]])
MD2 = np.array([[124,40],[124,30],[145,30],[145,40]])
MD = np.array([[125,38],[125,30],[145,30],[145,38]])
ST2 = np.array([[105,28],[105,18],[125,18],[125,28]])
#ST = np.array([[115,28],[115,22],[123,22],[123,28]])
ST = np.array([[113,28],[113,22],[123,22],[123,28]])
TEST=np.array([[121,35],[121,30],[125,30],[125,35]])

dcrds = {"NJ":NJ,"SJ":SJ, "KR":KR,"EC":EC,"SC":SC, "IC":IC, "PH":PH, "MD":MD, "ST":ST, "TEST":TEST}
#lname = ["NJ","SJ", "KR", "EC","SC", "IC", "PH"]
lname = ["MD","ST"]
#lname = ["TEST"]

lcrds = [dcrds[name] for name in lname]
rmasks = regionmask.Regions(lcrds, names=lname)

masks = rmasks.mask(a1loncnt, a1latcnt)
masks = ma.masked_invalid(masks)
print(masks.shape)

print(rmasks)
fig = plt.figure(figsize=(6,4))
axmap = fig.add_axes([0.1,0.1,0.8,0.8], projection=ccrs.PlateCarree())
axmap.set_extent([lllon,urlon,lllat,urlat])

#m = axmap.pcolormesh(a1lonbnd, a1latbnd, masks, transform=ccrs.PlateCarree())
m = axmap.contourf(a1loncnt, a1latcnt, masks)
rmasks.plot_regions(ax=axmap, add_label=False)

gl        = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,901, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
axmap.coastlines(color="k")

plt.show()

#-- Save masks ---
r3mask = rmasks.mask_3D(a1loncnt, a1latcnt)
for rname in lname:
    a2mask = r3mask.isel(region=(r3mask.names==rname)).values[0].astype("int16")

    maskdir = "/home/utsumi/temp/bams2020/mask"
    maskpath= maskdir + "/mask.aphro.%s.npy"%(rname)
    np.save(maskpath, a2mask)
    print(maskpath)


# %%

# %%
