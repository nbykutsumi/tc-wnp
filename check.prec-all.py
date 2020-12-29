import glob

print("")
print("*******************")
print("HPB")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/mqp-prec-all/HPB.*/prec.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)


print("")
print("*******************")
print("HPB_NAT")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/mqp-prec-all/HPB_NAT.*/prec.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)



