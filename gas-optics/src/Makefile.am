AM_CPPFLAGS = -I${top_srcdir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
if ENABLE_CUDA
cuda_objects = cuda_kernels.c
endif
lib_LIBRARIES = libgas_optics.a
libgas_optics_a_SOURCES = cfcs.c collision_induced_absorption.c gas_optics.c kernel_utils.c \
                          kernels.c launch.c molecules.c ozone_continuum.c \
                          parse_HITRAN_file.c RFM_voigt.c spectral_bin.c tips2017.c \
                          water_vapor_continuum.c $(cuda_objects)
if ENABLE_CUDA
gpu_object = gas_optics_gpu.$(OBJEXT)
objects = $(libgas_optics_a_OBJECTS)
libgas_optics_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif

include_HEADERS = cfcs.h collision_induced_absorption.h gas_optics.h molecules.h \
                  ozone_continuum.h spectral_bin.h water_vapor_continuum.h
