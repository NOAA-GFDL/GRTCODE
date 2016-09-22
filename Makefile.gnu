CC=gcc

CFLAGS=-Wall -Wextra -Wno-missing-field-initializers -g -O3 -fPIC -std=gnu99
CLIBS= -lnetcdf -lhdf5 -lhdf5_hl -lm
CFLAGS+= $(CLIBS)
#NETCDF_INCLUDES = -I${CPATH}
#INCLUDES = $(NETCDF_INCLUDES)
#CFLAGS += $(INCLUDES)
LDFLAGS= --shared


.PHONY: clean all

all: libgrtcode.so grtcode.x parseNetcdfRadiation.x TIPS_2011.unittest.x

libgrtcode.so: grtcode.o TIPS_2011.o parseHITRANfile.o parseNetcdfRadiation.o voigt.o outputNetcdfSpec.o continuum.o
	$(CC) $(LDFLAGS) $(CFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -DSKIPMAIN -c $< -o $@

%.o: %.cu
	$(CC) $(CFLAGS) -DSKIPMAIN -c $< -o $@

TIPS_2011.unittest.x: 
	$(CC) $(CFLAGS) TIPS_2011.c -o $@

grtcode.x: grtcode.o TIPS_2011.o parseHITRANfile.o parseNetcdfRadiation.o voigt.o outputNetcdfSpec.o continuum.o
	$(CC) $(CFLAGS)  $^ -o $@

parseNetcdfRadiation.x: parseNetcdfRadiation.c
	$(CC) $(CFLAGS) $^ -o $@
clean:
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
