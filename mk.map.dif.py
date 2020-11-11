from numpy import *
import myfunc.util as util
import sys, os
import numpy as np

ny,nx = 128, 256
lens = range(50)

prj     = "HAPPI"
model   = "MIROC5"
#lexpr    = ['ALL','P15','P20']
lexpr    = ['P20']
lens    = range(1,50)
#lens    = range(1,3)
nens    = len(lens)
res     = "128x256"
ny, nx  = 128, 256
miss    = -9999.
#thpr = 'p99.900'
thpr = '001.0'

dpintPath = {}

for expr in lexpr:
    if expr in ['P15','P20']:
        iYM  = [2106,1]
        eYM  = [2115,12]
    elif expr in ['ALL']:
        iYM  = [2006,1]
        eYM  = [2015,12]

    lYM  = util.ret_lYM(iYM,eYM)

    #** Initialize ********
    a3pint = np.zeros([ny,nx,nens],float32)
    #----------------------
    for iens,ens in enumerate(lens):
        a2sum = np.zeros([ny,nx],float32)
        a2num = np.zeros([ny,nx],int32)
        for [Year,Mon] in lYM:
            #srcDir = '/home/utsumi/mnt/lab_tank/utsumi/HAPPI/anlWS/tagpr/MIROC5/C20-P20-001/th.p99.900/210605'
            srcDir = '/home/utsumi/mnt/lab_tank/utsumi/HAPPI/anlWS/tagpr2/MIROC5/C20-%s-%03d/th.%s/%04d%02d'%(expr, ens, thpr, Year, Mon)
            #srcDir = '/media/disk2/share/HAPPI/anlWS/tagpr/MIROC5/C20-%s-%03d/th.%s/%04d%02d'%(expr, ens, thpr, Year, Mon)
            a2sumtmp = fromfile(srcDir + '/sum.tc.128x256', float32).reshape(ny,nx)
            a2numtmp = fromfile(srcDir + '/num.tc.128x256', float32).reshape(ny,nx)

            a2sum = a2sum + a2sumtmp
            a2num = a2num + a2numtmp

        a2pint = (ma.masked_where(a2num==0, a2sum)/a2num).filled(-9999.) * 60*60
        a3pint[iens] = a2pint

    #** Save data ********
    tmpDir = '/home/utsumi/temp/bams2020'
    util.mk_dir(tmpDir)

    pintPath = tmpDir + '/pintPath.%s.tc.npy'%(expr)
    dpintPath[expr] = pintPath

    np.save(pintPath, a3pint)
    print pintPath






