if USE_SINGLE_PRECISION
AM_CPPFLAGS = -DSINGLE_PRECISION
endif
lib_LIBRARIES = libargparse.a libgrtcode_utilities.a
libargparse_a_SOURCES = argparse.c
libgrtcode_utilities_a_SOURCES = curtis_godson.c device.c optics.c parse_csv.c \
                                 spectral_grid.c utilities.c verbosity.c
if ENABLE_CUDA
argparse.$(OBJEXT): argparse.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ -c $<
gpu_object = grtcode_utilities_gpu.$(OBJEXT)
objects = $(libgrtcode_utilities_a_OBJECTS)
libgrtcode_utilities_a_LIBADD = $(gpu_object)
include ${top_srcdir}/cuda.mk
endif

include_HEADERS = device.h extern.h floating_point_type.h grtcode_utilities.h optics.h \
                  spectral_grid.h utilities.h verbosity.h
