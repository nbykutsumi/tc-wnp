import glob

print("")
print("*******************")
print("HPB")
lsrcPath = sorted(glob.glob("/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM/XX-HPB-*/6hr/cwise/2010/1230/dura*"))
for srcPath in lsrcPath:
    print(srcPath)

#
print("")
print("*******************")
print("HPB_NAT")
lsrcPath = sorted(glob.glob("/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM/XX-HPB_NAT-*/6hr/cwise/2010/1230/dura*"))
for srcPath in lsrcPath:
    print(srcPath)



