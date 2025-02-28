AC_PREREQ([2.59])
AC_INIT([grtcode], [0.0], [gfdl.climate.model.info@noaa.gov])
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign dist-zip subdir-objects])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_CC
AC_C_CONST
AC_PROG_CPP
AC_PROG_CXX
AM_PROG_AR
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_SYS_LARGEFILE

#Check for math library.
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR([cannot find math.h])])
AC_CHECK_LIB([m], [pow], [], [AC_MSG_ERROR([cannot link -lm.])])

#Check for netcdf.
AC_CHECK_HEADER([netcdf.h], [], [AC_MSG_ERROR([cannot find netcdf.h])])
AC_CHECK_LIB([netcdf], [nc_open], [], [AC_MSG_ERROR([could not link libnetcdf.])])

#Provide an option to build in single precision.
AC_ARG_ENABLE([single_precision],
[--enable-single-precision Use single precision],
[case "${enableval}" in
  yes) single_precision=true ;;
  no)  single_precision=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-single-precision]) ;;
esac],[single_precision=false])
AM_CONDITIONAL([USE_SINGLE_PRECISION], [test x$single_precision = xtrue])

#Provide an option to build with CUDA.
AC_ARG_VAR([NVCC], [nvcc compiler command])
AC_ARG_VAR([NVCCFLAGS], [nvcc compiler flags])
AC_ARG_VAR([NVCCLDFLAGS], [nvcc linker flags])
AC_ARG_ENABLE([cuda],
[--enable-cuda Enable ability to run on NVIDIA GPUs.],
[case "${enableval}" in
  yes) cuda=true ;;
  no)  cuda=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-cuda]) ;;
esac],[cuda=false])
AM_CONDITIONAL([ENABLE_CUDA], [test x$cuda = xtrue])
if test x$cuda = xtrue
then
  AC_PATH_PROG([NVCC], [nvcc], [no])
  if test x$NVCC != xno
  then
    AC_MSG_CHECKING([whether nvcc works])
    cat > conftest.cu <<EOF
__global__ void test(void)
{
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    return;
}
int main(void)
{
    test<<<1, 32, 0, 0>>>();
    return 0;
}
EOF
    if $NVCC -x cu -w -o conftest.o -c conftest.cu && \
       $NVCC -dlink -o conftest-device.o conftest.o && \
       $CXX $LDFLAGS -o conftest.x conftest.o conftest-device.o -lcudart && \
       ./conftest.x
    then
      AC_MSG_RESULT([yes])
      rm -f conftest.cu conftest.x
    else
      AC_MSG_RESULT([no])
      rm -f conftest.cu conftest.x
      AC_MSG_ERROR([nvcc could not compile/run executable])
    fi
  else
    AC_MSG_ERROR([could not find nvcc])
  fi
fi

#Provide an option to build fortran bindings.
AC_ARG_ENABLE([fortran],
[--enable-fortran Build fortran bindings.],
[case "${enableval}" in
  yes) fortran=true ;;
  no)  fortran=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-fortran]) ;;
esac],[fortran=false])
AM_CONDITIONAL([FORTRAN_BINDINGS], [test x$fortran = xtrue])
AC_PROG_FC

#Provide an option to use DISORT.
AC_ARG_ENABLE([disort],
[--enable-disort Use disort solver.],
[case "${enableval}" in
  yes) disort=true ;;
  no)  disort=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-disort]) ;;
esac],[disort=false])
AM_CONDITIONAL([USE_DISORT], [test x$disort = xtrue])
if test x$disort = xtrue
then
  AC_CHECK_HEADER([cdisort.h], [], [AC_MSG_ERROR([cannot find cdisort.h])])
  AC_CHECK_LIB([cdisort], [c_disort], [], [AC_MSG_ERROR([could not link libcdisort.])], [-lm])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 utilities/Makefile utilities/src/Makefile
                 gas-optics/Makefile gas-optics/src/Makefile
                 shortwave/Makefile shortwave/src/Makefile shortwave/test/Makefile
                 longwave/Makefile longwave/src/Makefile longwave/test/Makefile
                 fortran-bindings/Makefile])
AC_OUTPUT
