AM_CPPFLAGS = -I${top_srcdir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
lib_LIBRARIES = liblongwave.a
liblongwave_a_SOURCES = longwave.c
if ENABLE_CUDA
gpu_object = longwave_gpu.$(OBJEXT)
objects = $(liblongwave_a_OBJECTS)
liblongwave_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif

include_HEADERS = longwave.h
