AM_CPPFLAGS = -I${top_srcdir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
if USE_DISORT
AM_CPPFLAGS += -DUSE_DISORT
endif
lib_LIBRARIES = libshortwave.a
libshortwave_a_SOURCES = disort_shortwave.c rayleigh.c shortwave.c solar_flux.c
if ENABLE_CUDA
gpu_object = shortwave_gpu.$(OBJEXT)
objects = $(libshortwave_a_OBJECTS)
libshortwave_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif

include_HEADERS = disort_shortwave.h rayleigh.h shortwave.h solar_flux.h
