#!/bin/bash -e

source $MODULESHOME/init/bash

module unload cray-netcdf cray-hdf5 fre
module unload PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module load PrgEnv-intel/8.6.0
module unload intel intel-classic intel-oneapi
module load intel-classic/2023.2.0
module load cray-hdf5/1.12.2.11
module load cray-netcdf/4.9.0.9
module unload cray-libsci
module load cray-libsci/24.11.0

export KMP_STACKSIZE="512m"
export NC_BLKSZ="1M"
export F_UFMTENDIAN="big"

PREFIX="/gpfs/f5/gfdl_m/scratch/Jing.Feng/line-by-line"
CC="cc"
CPPFLAGS="-I${PREFIX}/include"
CFLAGS="-g -O3 -qopenmp"
MPICC="cc"
CXX="CC"
CXXFLAGS="-g -O3 -qopenmp"
MPICXX="CC"
LDFLAGS="-L${PREFIX}/lib"

#debug
#CPPFLAGS="-I${PREFIX}/include"
#CFLAGS="-g -O0 -fno-omit-frame-pointer -qopenmp"
#CXXFLAGS="-g -O0 -fno-omit-frame-pointer -qopenmp"
