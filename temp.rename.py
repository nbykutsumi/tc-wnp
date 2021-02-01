import os, shutil

radkm = 500
lscen = ["HPB","HPB_NAT"]
iY,eY = 1990,2010
lens  = range(1,50+1)

for scen in lscen:
    for ens in lens:
        precbasedir = "/home/utsumi/temp/bams2020/tc-prec-%04dkm"%(radkm)
        avedir = precbasedir + "/ens-ave-%04d-%04d"%(iY,eY)
        #ipath= avedir + "/prec-ave.%s.%03d.npy"%(scen,ens)
        #opath= avedir + "/prec-tc-ave.%s.%03d.npy"%(scen,ens)
        #print(ipath)
        #print(os.path.exists(ipath))
        #os.rename(ipath,opath)
