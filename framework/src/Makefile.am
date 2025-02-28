AM_CPPFLAGS = -I${top_srcdir}/shortwave/src \
              -I${top_srcdir}/longwave/src \
              -I${top_srcdir}/gas-optics/src \
              -I${top_srcdir}/utilities/src \
              -I${top_srcdir}/dependencies/clouds-lib/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif
noinst_LIBRARIES = libdriver.a
libdriver_a_SOURCES = driver.c
