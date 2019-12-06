AM_CPPFLAGS = -I${top_builddir}/utilities/src \
              -I${top_builddir}/gas-optics/src \
              -I${top_builddir}/shortwave/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
if FORTRAN_BINDINGS
lib_LIBRARIES = libgrtcode_fortran.a
libgrtcode_fortran_a_SOURCES = malloc_structs.c grtcode_fortran.F90
grtcode.mod: grtcode_fortran.$(OBJEXT)
CLEANFILES = grtcode.mod
include_HEADERS = grtcode.mod
if ENABLE_CUDA
gpu_object = malloc_structs_gpu.$(OBJEXT)
objects = malloc_structs.$(OBJEXT)
libgrtcode_fortran_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif
endif
