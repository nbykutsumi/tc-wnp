import os, sys, shutil, glob

lscen = ['HPB','HPB_NAT']
#lscen = ['HPB']
lens = range(1,10+1)
#lens = [2]

#for scen in lscen:
#    for ens in lens:
#        #ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-*/6hr/*/*/*/*320x640'
#        ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-%s-%03d/6hr/*/*/*/*320x640'%(scen,ens)
#        lpath = glob.glob(ssearch)
#       
#        for spath in lpath: 
#            print spath
#            os.remove(spath)    

#for scen in lscen:
#    for ens in lens:
#        #ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-*/6hr/*/*/*/*320x640'
#        ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-%s-%03d/6hr/clist/*/*/*bn'%(scen,ens)
#        lpath = glob.glob(ssearch)
#        print len(lpath) 
#        for spath in lpath: 
#            print spath
#            os.remove(spath)    
#

#for scen in lscen:
#    for ens in lens:
#        #ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-*/6hr/*/*/*/*320x640'
#        ssearch = '/tank/utsumi/WS/d4PDF_GCM/XX-%s-%03d/6hr/pgrad'%(scen,ens)
#        lpath = glob.glob(ssearch)
#        print len(lpath)
#        for spath in lpath: 
#            print spath
#            shutil.rmtree(spath)



