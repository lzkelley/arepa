import numpy as np

Base="./output/"
SnapSkipFac=1
num=40

fname=Base + "/treedata/sub_desc_sf"+str(SnapSkipFac)+"_"+str(num).zfill(3)

f=open(fname, "r")
TotNsubhalos = np.fromfile(f, dtype=np.uint32, count=1)[0]
descendant_halonr = np.fromfile(f, dtype=np.uint32, count=TotNsubhalos)
descendant_snapnr = np.fromfile(f, dtype=np.uint32, count=TotNsubhalos)
f.close()


print TotNsubhalos
print descendant_halonr.shape
print descendant_snapnr.shape

print descendant_halonr
print descendant_snapnr
