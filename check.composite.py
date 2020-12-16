import glob

print("")
print("---------------------------------------")
print("HPB")
print("---------------------------------------")
ssearch = "/home/utsumi/temp/bams2020/composite/sst--9999.ex-3.00.tc-3.00.wc-0.0-wind--9999-wdif--9999-du-36/HPB.*/*/dura.*.12.npy"
lsrcpath = sorted(glob.glob(ssearch))
for srcpath in lsrcpath:
    print(srcpath)

print("")
print("---------------------------------------")
print("HPB_NAT")
print("---------------------------------------")
ssearch = "/home/utsumi/temp/bams2020/composite/sst--9999.ex-3.00.tc-3.00.wc-0.0-wind--9999-wdif--9999-du-36/HPB_NAT.*/*/dura.*.12.npy"
lsrcpath = sorted(glob.glob(ssearch))
for srcpath in lsrcpath:
    print(srcpath)



print("")
print("---------------------------------------")
print("OBS")
print("---------------------------------------")
ssearch = "/home/utsumi/temp/bams2020/composite/obs-era5-tc/*/prec000.*.12.npy"
lsrcpath = sorted(glob.glob(ssearch))
for srcpath in lsrcpath:
    print(srcpath)
