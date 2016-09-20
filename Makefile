CC=cc
NVCC=nvcc

CFLAGS=' -Wall -Wextra -Wno-missing-field-initializers -g -O3 '
CLIBS= -lnetcdf -lhdf5 -lhdf5_hl
CFLAGS+= $(CLIBS)
#NETCDF_INCLUDES = -I${CPATH}
#INCLUDES = $(NETCDF_INCLUDES)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
#LDFLAGS= --shared
NVCCFLAGS= --compiler-options $(CFLAGS) -g -O3 --use_fast_math 
NVCCLDFLAGS= --compiler-options '-shared -fPIC '
NVOPT= --maxrregcount=32
NVCCFLAGS += $(NVOPT)



.PHONY: clean all

all: libparseNetcdfRadiation.so libgrtcode.so grtcode.x parseNetcdfRadiation.x TIPS_2011.unittest.x

libgrtcode.so: grtcode.o TIPS_2011.o parseHITRANfile.o cudaHelpers.o parseNetcdfRadiation.o voigt.o outputNetcdfSpec.o continuum.o
	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

libparseNetcdfRadiation.so: parseNetcdfRadiation.o
	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

%.o: %.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

%.o: %.cu
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

TIPS_2011.unittest.x: 
	$(NVCC) -x cu $(NVCCFLAGS) TIPS_2011.c cudaHelpers.cu -o $@

grtcode.x: grtcode.o TIPS_2011.o parseHITRANfile.o cudaHelpers.o parseNetcdfRadiation.o voigt.o outputNetcdfSpec.o continuum.o
	$(NVCC) $(NVCCFLAGS)  $^ -o $@

parseNetcdfRadiation.x: parseNetcdfRadiation.c
	$(NVCC) $(NVCCFLAGS) $^ -o $@
clean:
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
