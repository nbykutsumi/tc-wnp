import glob

print("")
print("*******************")
print("HPB")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/tc-prec-annmax/HPB.*/num.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)


print("")
print("*******************")
print("HPB_NAT")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/tc-prec-annmax/HPB_NAT.*/num.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)



