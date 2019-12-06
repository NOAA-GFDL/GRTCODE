AM_CPPFLAGS = -I${top_builddir}/utilities/src
if USE_SINGLE_PRECISION
AM_CPPFLAGS += -DSINGLE_PRECISION
endif

check_PROGRAMS = check_results
check_results_SOURCES = check_results.c
check_results_LDADD = ${top_builddir}/utilities/src/libgrtcode_utilities.a \
                      ${top_builddir}/utilities/src/libargparse.a \
                      -lnetcdf -lm
if ENABLE_CUDA
check_results_LDADD += -lcudart
endif
TESTS = test_circ
