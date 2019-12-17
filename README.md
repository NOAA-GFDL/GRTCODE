[![Build Status](https://travis-ci.com/NOAA-GFDL/GRTCODE.svg?branch=master)](https://travis-ci.com/NOAA-GFDL/GRTCODE)

# GPU-able Radiative Transfer code (GRTCODE)
Fork of [GRTCODE](https://gitlab.com/geebdubya/GRTCODE)
originally written by Garret Wright ([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3451457.svg)](https://doi.org/10.5281/zenodo.3451457)).


# Overview
GRTCODE consists of a column-based line-by-line radiative transfer
model capable of running on NVIDIA GPUs via CUDA-c and many-core CPUs
via OpenMP.  This package consists of several general purpose libraries
and a set of binaries that specifically target popular radiative transfer
benchmarks.


### Libraries

##### libgrtcode_utilities.a
A utility library containing helper functions, error handlers, etc.

##### libgas_optics.a
A library that calculates the optical properties (i.e. optical depths) of an
atmospheric column assumed to be made up of a series of "pressure levels".
This library uses the HITRAN database to calculate spectra for gases including
water vapor, ozone, carbon dioxide, etc. and several species of CFCs, HCFCs, and HFCs.

##### liblongwave.a
A four-stream longwave solver that uses a 2-term Pade approximation.

##### libshortwave.a
A two-stream solver that uses the delta-Eddington and Adding methods.

##### libgrtcode_fortran.a
FORTRAN bindings to the above libraries.

### Binaries

##### circ
Program that calculates longwave and shortwave radiative fluxes for the
[CIRC cases](https://circ.gsfc.nasa.gov/index.html).

##### rfmip-irf
Program that calculates longwave and shortwave radiative fluxes for the
[RFMIP-IRF offline benchmark cases](https://doi.org/10.5194/gmd-9-3447-2016).

# Requirements
This repository is setup to use the [GNU Build System](https://www.gnu.org/software/automake/manual/html_node/GNU-Build-System.html).
The libraries and binaries are written in c (ISO/IEC 9899:1999 standard) and thus
requires a c compiler, such as the freely available gcc.  In order to run on NVIDIA GPUs,
a c++ compiler (such as g++), CUDA, and the NVCC compiler are also required.
The optional FORTRAN front-end library obviously also requires a FORTRAN
compiler (such as gfortran).  Some of the binaries also require the netCDF library.
A tarball with example input data can be downloaded from
[GFDL's FTP server](ftp://ftp2.gfdl.noaa.gov/perm/GFDL_pubrelease/test_data/grtcode-data.tar.gz).
To run the tests, download the example data, untar in the base of the repository,
and then run ```make check``` after building.

# Building

### CPU-only
To build using default settings, run the normal GNU Build System commands:

```
$ autoreconf --install
$ ./configure [--prefix <where_to_install>] [--enable-single-precision]
$ make
$ make check (this is optional but recommended)
$ make install
```

As usual, the default compilers and flags can be
overridden by specifying CC, CPPFLAGS, CFLAGS, LDFLAGS) when running
configure in the usual fashion:

```
$ ./configure CC=gcc CFLAGS='-O3 -fopenmp' LDFLAGS='-fopenmp'
```

This library can take advantage of thread-level parallelism through the
use of OpenMP (the default settings run single-threaded).  If your compiler
supports OpenMP, make sure you activate it by setting the corresponding
compiler option.

### With CUDA
To build a CUDA-enabled version, run:

```
$ autoreconf --install
$ ./configure --enable-cuda
$ make
$ make check (this is optional, but recommended and requires downloading/untar-ing the example input data.)
$ make install
```

As stated above, CUDA-enable builds require NVCC and a c++ compiler.  If not
included in your $PATH, these (as well as NVCCFLAGS and NVCCLDFLAGS) can be
set when running configure:

```
$ ./configure --enable-cuda CXX=g++ NVCC=/usr/local/cuda/bin/nvcc NVCCFLAGS='-g -O3' \
              NVCCLDFLAGS=-L/usr/local/cuda/lib
```

### FORTRAN bindings
A FORTRAN front-end can be built for either the CPU-only or CUDA-enabled versions
by running:

```
$ autoreconf --install
$ ./configure --enable-fortran
$ make
$ make check (this is optional, but recommended and requires downloading/untar-ing the example input data.)
$ make install
```

The default FORTRAN compiler and flags can also be overridden in the usual fashion:

```
./configure --enable-fortran FC="gfortran" FCFLAGS="-g -O3"
```
