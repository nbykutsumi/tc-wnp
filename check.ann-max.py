import glob

print("")
print("*******************")
print("HPB")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/max-prec-d/HPB.*/ann-max-prec.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)


print("")
print("*******************")
print("HPB_NAT")
lsrcPath = sorted(glob.glob("/home/utsumi/temp/bams2020/max-prec-d/HPB_NAT.*/ann-max-prec.*.npy"))
for srcPath in lsrcPath:
    print(srcPath)



