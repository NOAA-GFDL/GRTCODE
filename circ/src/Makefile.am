AM_CPPFLAGS = -I${top_srcdir}/framework/src \
              -I${top_srcdir}/shortwave/src \
              -I${top_srcdir}/longwave/src \
              -I${top_srcdir}/cloud-optics/src \
              -I${top_srcdir}/gas-optics/src \
              -I${top_srcdir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
AM_LDFLAGS = -L${top_builddir}/dependencies/clouds-lib/.libs
bin_PROGRAMS = circ
circ_SOURCES = circ.c
circ_LDADD = ${top_builddir}/framework/src/libdriver.a \
             ${top_builddir}/shortwave/src/libshortwave.a \
             ${top_builddir}/longwave/src/liblongwave.a \
             ${top_builddir}/gas-optics/src/libgas_optics.a \
             ${top_builddir}/utilities/src/libgrtcode_utilities.a \
             ${top_builddir}/utilities/src/libargparse.a \
             -lcloudsc -lnetcdf
if ENABLE_CUDA
circ_LDADD += -lcudart
nodist_EXTRA_circ_SOURCES = dummy.cpp
endif
