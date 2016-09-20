#!/usr/bin/env python

##This is just a hack to take certain collections of flat files and write them in a single netcdf file.
### Note you will want to change NETCDF4 to something else, as that format breaks some tools, if you forget use nccopy to an classic 64bit offst.

import sys
import numpy as np
import glob
#ugh
from netCDF4 import Dataset


def getFnames(path,ext='spc'):
    return glob.glob(path+'/*-*.'+ext)

def parsefn(fn,path=None):
    pre = 0
    if path is not None:
        pre = len(path)
    lat,lon = fn[pre:].split('.')[0].split('-')
    return int(lat),int(lon)

def parseDavidSpcFile(fn,nhdr=1):    
    with open(fn,'r') as f:
        lines = f.readlines()        
    hdr = lines[0]
    skip = 1 #1 if you start at 0wavenumber..., else 0
    lines = lines[1+skip:]
    n = len(lines)
    nleadingCol=1
    ntrailingCol=0#1    
    nc = len(lines[100].strip().split()[nleadingCol:]) - ntrailingCol
    V = np.empty(n,dtype=np.float64)
    D = np.empty((nc,n),dtype=np.float64)
    for l in xrange(n):
        line = lines[l].strip().split()
        V[l] = line[0]
        for col in xrange(nc):
            #nc-col-1 reverse z order
            D[nc-col-1][l] = line[1+col]

    return hdr,V,D

                                                                                    
if __name__ == "__main__":

    NLAT=90
    NLON=144
    NZ=48
    NWV=3000
    
    assert len(sys.argv) == 3
        
    outF = sys.argv[1]
    inPath = sys.argv[2]

    rootgrp = Dataset(outF, "w", format="NETCDF4")
    level = rootgrp.createDimension("pfull", NZ)
    time = rootgrp.createDimension("time", None)
    lat = rootgrp.createDimension("lat", NLAT)
    lon = rootgrp.createDimension("lon", NLON)
    wv = rootgrp.createDimension("wavenumber", NWV)

    opt = rootgrp.createVariable("OpticalDepth", "f4", ("time", "lat", "lon", "pfull", "wavenumber",))
    
    names = getFnames(inPath)
    for nm in names[:]:
        lat,lon = parsefn(nm, inPath)
        print lat, lon
        hdr, v, d = parseDavidSpcFile(nm)
        print d.shape
        assert d.shape == (NZ,NWV)
        opt[0, lat, lon, :, :] = d


    

    rootgrp.close()
