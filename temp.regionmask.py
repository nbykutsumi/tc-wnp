# %%

import numpy as np
from numpy import ma
import regionmask
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib
import d4PDF
matplotlib.use("Agg")
%matplotlib inline

a1latcnt = d4PDF.Lat()
a1loncnt = d4PDF.Lon()
a1latbnd = d4PDF.LatBnd()
a1lonbnd = d4PDF.LonBnd()

[[lllat,lllon],[urlat,urlon]] = [[0,100],[50,180]]   

JP = np.array([[123,25],[123,23],[130,23],[135,30],[142,30],[150,46],[140,46],[135,40],[129,34]])   # [lon, lat]

lcrds = [JP]
lname = ["JP"]
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

gl        = axmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, linestyle=":", color="k", alpha=0.8)
xticks   = np.arange(-180, 180+1, 15)
yticks   = np.arange(-90,901, 15)
gl.xlocator = mticker.FixedLocator(xticks)
gl.ylocator = mticker.FixedLocator(yticks)
axmap.coastlines(color="k")

plt.show()
# %%
