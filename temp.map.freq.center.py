# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
#----------------------------------
import sys, os, pickle
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   detect_fsub import *
from   datetime import datetime, timedelta
from   importlib import import_module
import numpy as np
import matplotlib.pyplot as plt
import util
import d4PDF
import scipy.stats as stats
#--------------------------------------
calcobj= True
#calcobj= False
figobj = True
#figobj = False

iDTime = datetime(2000,1,1,0)
eDTime = datetime(2000,12,31,18)
#eDTime = datetime(2010,12,31,18)
#eDTime = datetime(2000,1,31,0)
dDTime = timedelta(hours=6)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
#-----------------
prj     = "d4PDF"
model   = "__"
expr    = 'XX'
scen    = 'HPB'
lens    = [1]
res     = "320x640"
noleap  = False
vname   = 'pgrad'
dattype = {'pgrad':'float32', 'ipos':'int32'}[vname]
wsbasedir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

detectName = 'wsd_d4pdf_20200428'
figdir  = '/home/utsumi/temp/bams2020'
#----------------------------------
a1latin = d4PDF.Lat()
a1lonin = d4PDF.Lon()
nyin    = len(a1latin)
nxin    = len(a1lonin)

miss_int= -9999

BBox  = [[0,100],[45,180]]
dgrid = 5
a1latbnd = np.arange(-90, 90+0.01, dgrid)
a1lonbnd = np.arange(0, 360+0.01, dgrid)
ny = a1latbnd.shape[0] -1
nx = a1lonbnd.shape[0] -1
#**************************
def draw_map(a2dat, dpara):
    figmap   = plt.figure(figsize=(6,4))
    axmap    = figmap.add_axes([0.1, 0.1, 0.7, 0.8])
    [[lllat,lllon],[urlat,urlon]] = dpara['BBox']
    M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

    #-- pcolormesh --------------
    mycm = 'gist_stern_r'
    vmin, vmax = 0, 0.02
    #vmin, vmax = 0, 1.0
    X,Y = np.meshgrid(a1lonbnd,a1latbnd)
    im  = M.pcolormesh(X,Y,a2dat, cmap=mycm, vmin=vmin,vmax=vmax)
    #-- Colorbar ------
    cax = figmap.add_axes([0.82,0.2,0.05,0.6])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)

    #-- coastline ---------------
    print "coastlines"
    M.drawcoastlines(color="k")

    #-- Merdians and Parallels --
    M.drawmeridians(np.arange(-180,180+1,15), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50, yoffset=2)
    M.drawparallels(np.arange(-60,60+1,15), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    #-- Tiele -----
    stitle = dpara['title']
    axmap.set_title(stitle)
    #-- Save ------
    figpath = dpara['figpath']
    plt.savefig(figpath)
    print figpath



#************************************
# d4PDF
#************************************
for ens in lens:
    if calcobj is True:
        a2dat = np.zeros([320,640], 'int32')
        for DTime in lDTime:
            print ens,DTime
            Year,Mon,Day,Hour = DTime.timetuple()[:4]
            srcdir = wsbasedir + '/%s-%s-%03d/6hr/%s/%04d/%02d'%(expr,scen,ens,vname,Year,Mon)
            srcpath = srcdir + '/%s.%04d%02d%02d%02d.320x640'%(vname,Year,Mon,Day,Hour)
            a2in = np.fromfile(srcpath, 'float32').reshape(320,640)
            a2in = ma.masked_greater(a2in,0).mask.astype('int32')
            a2dat = a2dat + a2in


        #--- Regrid ---                 
        a2lonin, a2latin = np.meshgrid(a1lonin,a1latin) 
        statistic, x_edge, y_edge, binnumber = stats.binned_statistic_2d(x=a2latin.flatten(), y=a2lonin.flatten(), values=a2dat.flatten(), statistic='sum', bins=[a1latbnd, a1lonbnd])

        a2out = statistic
        a2out = a2out.astype('float32')/len(lDTime)

    #-- Save --------
    savedir = figdir + '/pickle'
    util.mk_dir(savedir)
    outpath = savedir + '/map.freq.%s.%s.%s-%s-%03d.npy'%(prj,vname,expr,scen,ens)
    if calcobj is True:
        np.save(outpath, a2out)

    #--- Figure ------
    if figobj is True:
        a2fig   = np.load(outpath)
        run = '%s-%s-%03d'%(expr,scen,ens)
        figpath = figdir + '/map.freq.%s.%s.%s-%s-%03d.png'%(prj,vname,expr,scen,ens)
        dpara = {}
        dpara['title'] = run
        dpara['figpath'] = figpath
        dpara['BBox'] = BBox
        draw_map(a2fig, dpara)




# %%
