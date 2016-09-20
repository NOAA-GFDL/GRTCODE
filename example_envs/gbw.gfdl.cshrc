set echo

module load fre/bronx-10

#CUDA
module unload cuda
setenv GPU "/net/gbw/gpgpu"
setenv PATH "/net/gbw/gpgpu/cuda/bin:${PATH}"
setenv LD_LIBRARY_PATH "/net/gbw/gpgpu/cuda/lib64:${LD_LIBRARY_PATH}" 

unset echo

module list

