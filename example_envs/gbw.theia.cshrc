set echo

module unload impi openmpi intel pgi
module load cuda/8.0
setenv LD_LIBRARY_PATH /scratch4/GFDL/gfdlscr/Garrett.Wright/local/netcdf-4.4.1/lib/:/scratch4/GFDL/gfdlscr/Garrett.Wright/local/hdf5-1.8.9/lib/:$LD_LIBRARY_PATH

unset echo
module list

