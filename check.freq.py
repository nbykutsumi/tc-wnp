import glob

print("")
print("*******************")
print("HPB")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/map-freq/sst-270.ex-3.00.tc-3.00.wc-0.0-wind-12-wdif--9999-du-36/HPB-*/a2count.tc.obj.2010.npy"))
for srcPath in lsrcPath:
    print(srcPath)


print("")
print("*******************")
print("HPB_NAT")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/map-freq/sst-270.ex-3.00.tc-3.00.wc-0.0-wind-12-wdif--9999-du-36/HPB_NAT-*/a2count.tc.obj.2010.npy"))
for srcPath in lsrcPath:
    print(srcPath)



