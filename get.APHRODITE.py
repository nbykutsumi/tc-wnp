import os, sys
import subprocess

hostname = "http://aphrodite.st.hirosaki-u.ac.jp"
uid   = "nobuyuki2549@gmail.com"
pswd = "SGLbu6eB"

##-- V1101 ---
#srcdir = "/product/APHRO_V1101/APHRO_MA/050deg"
#outdir = "/tank/utsumi/data/APHRO/APHRO_V1101/APHRO_MA/050deg"
#iYear = 1951
#eYear = 2007
##eYear = 1951
#lYear = list(range(iYear,eYear+1))
#for Year in lYear:
#    fname = "APHRO_MA_050deg_V1101.%04d.gz"%(Year)
#    scmd  = 'wget --http-user=%s --http-password=%s -r --no-directories --no-parent -A %s %s/%s/ --directory-prefix=%s'%(uid, pswd, fname, hostname, srcdir, outdir)
#    lcmd  = scmd.split()
#    subprocess.call(lcmd)
#
##-- ctrl file --
#fname = "APHRO_MA_050deg_V1101.ctl.gz"
#scmd  = 'wget --http-user=%s --http-password=%s -r --no-directories --no-parent -A %s %s/%s/ --directory-prefix=%s'%(uid, pswd, fname, hostname, srcdir, outdir)
#lcmd  = scmd.split()
#subprocess.call(lcmd)


#-- V1101EX_R1 ---
srcdir = "/product/APHRO_V1101EX_R1/APHRO_MA/050deg"
outdir = "/tank/utsumi/data/APHRO/APHRO_V1101EX_R1/APHRO_MA/050deg"
iYear = 2008
eYear = 2015
#eYear = 1951
lYear = list(range(iYear,eYear+1))
for Year in lYear:
    fname = "APHRO_MA_050deg_V1101_EXR1.%04d.gz"%(Year)
    scmd  = 'wget --http-user=%s --http-password=%s -r --no-directories --no-parent -A %s %s/%s/ --directory-prefix=%s'%(uid, pswd, fname, hostname, srcdir, outdir)
    lcmd  = scmd.split()
    subprocess.call(lcmd)

#-- ctrl file --
fname = "APHRO_MA_050deg_V1101_EXR1.ctl.gz"
scmd  = 'wget --http-user=%s --http-password=%s -r --no-directories --no-parent -A %s %s/%s/ --directory-prefix=%s'%(uid, pswd, fname, hostname, srcdir, outdir)
lcmd  = scmd.split()
subprocess.call(lcmd)



