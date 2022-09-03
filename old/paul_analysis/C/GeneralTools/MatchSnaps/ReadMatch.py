import readsubfHDF5
import numpy as np

Base1='/n/hernquistfs1/Illustris/Runs/Illustris-3/output/'
Base2='/n/hernquistfs1/Illustris/Runs/Illustris-Dark-3/output/'

MatchBase="./output/"
snapnum=100

fname=MatchBase + "/matchdata/sub_match_"+str(snapnum).zfill(3)


cat1 = readsubfHDF5.subfind_catalog(Base1, snapnum, keysel=["SubhaloPos","SubhaloLenType"])
cat2 = readsubfHDF5.subfind_catalog(Base2, snapnum, keysel=["SubhaloPos","SubhaloLenType"])


ch=200

print cat1.SubhaloPos[ch,:]
print cat2.SubhaloPos[match_halonr[ch],:]
