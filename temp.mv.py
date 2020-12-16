import shutil, os
import glob
import util
#lvname = ["dura","idate","ipos","land","lat","lon","prec","sp","time","tsum","vort","wind","wlw","wup"]
lens=range(1,20+1)
#lens=

basedir = "/home/utsumi/mnt/lab_tank/utsumi/hometemp/bams2020/composite/sst--9999.ex-3.00.tc-3.00.wc--9999.0-wind--9999-wdif--9999-du-36"
for ens in lens:
    ssearch = basedir + "/HPB.%03d/*.npy"%(ens)
    print(ssearch)
    lsrcpath = sorted(glob.glob(ssearch))
    print(ens, len(lsrcpath))
    for srcpath in lsrcpath:
        Year   = int(os.path.basename(srcpath).split(".")[1])
        srcdir = os.path.dirname(srcpath)
        outdir = srcdir + "/%04d"%(Year)
#        print(outdir)
#        util.mk_dir(outdir)
#        shutil.move(srcpath, outdir + "/")
