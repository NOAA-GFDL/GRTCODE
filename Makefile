ifndef CC
CC = gcc
endif
CPPFLAGS += -Iutilities/src -Igas-optics/src -Ishortwave/src -Ilongwave/src -Iclouds -Iframework/src -Itesting_harness/src
ifndef CFLAGS
CFLAGS = -g -O0 -Wall -Wextra -pedantic -fopenmp
endif

OBJECTS = build/longwave.o \
          build/disort_shortwave.o \
          build/rayleigh.o \
          build/shortwave.o \
          build/solar_flux.o \
          build/RFM_voigt.o \
          build/cfcs.o \
          build/collision_induced_absorption.o \
          build/cuda_kernels.o \
          build/gas_optics.o \
          build/kernel_utils.o \
          build/kernels.o \
          build/launch.o \
          build/molecules.o \
          build/ozone_continuum.o \
          build/parse_HITRAN_file.o \
          build/spectral_bin.o \
          build/tips2017.o \
          build/water_vapor_continuum.o \
          build/argparse.o \
          build/curtis_godson.o \
          build/device.o \
          build/optics.o \
          build/parse_csv.o \
          build/spectral_grid.o \
          build/utilities.o \
          build/verbosity.o

CLOUD_OBJECTS = build/clouds_lib.o \
                build/hu_stamnes.o \
                build/ice_cloud_optics.o \
                build/incomplete_beta.o \
                build/netcdf_utils.o \
                build/optics_utils.o \
                build/stochastic_cloud.o

TEST_OBJECTS = build/test_harness.o \
               build/test_curtis_godson.o \
               build/test_device.o \
               build/test_spectral_grid.o \
               build/test_verbosity.o \
               build/test_parse_csv.o \
               build/test_utilities.o \
               build/test_optics.o \
               build/test_longwave.o \
               build/test_shortwave.o \
               build/test_RFM_voigt.o \
               build/test_cfcs.o \
               build/test_collision_induced_absorption.o \
               build/test_gas_optics.o \
               build/test_kernel_utils.o \
               build/test_molecules.o \
               build/test_ozone_continuum.o \
               build/test_tips2017.o \
               build/test_water_vapor_continuum.o \
               build/test_parse_HITRAN_file.o \
               build/test_kernels.o \
               build/test_spectral_bin.o
#              build/test_disort_shortwave.o \
#              build/test_rayleigh.o \
#              build/test_solar_flux.o \
#              build/test_cuda_kernels.o \
#              build/test_launch.o \
#              build/test_argparse.o \

TESTS = build/test_curtis_godson \
        build/test_device \
        build/test_spectral_grid \
        build/test_verbosity \
        build/test_parse_csv \
        build/test_utilities \
        build/test_optics \
        build/test_longwave \
        build/test_shortwave \
        build/test_kernel_utils \
        build/test_molecules \
        build/test_cfcs \
        build/test_collision_induced_absorption \
        build/test_gas_optics \
        build/test_ozone_continuum \
        build/test_tips2017 \
        build/test_water_vapor_continuum \
        build/test_RFM_voigt \
        build/test_parse_HITRAN_file \
        build/test_kernels \
        build/test_spectral_bin
#       build/test_disort_shortwave \
#       build/test_rayleigh \
#       build/test_solar_flux \
#       build/test_cuda_kernels \
#       build/test_launch \
#       build/test_argparse \

all: build/libgrtcode.a

check: $(TESTS)
	-build/test_curtis_godson
	build/test_device
	build/test_spectral_grid
	build/test_parse_csv

circ: build/circ

era5: build/era5

rfmip-irf: build/rfmip-irf

# Core libraries.
build/libgrtcode.a: $(OBJECTS)
	ar rcs $@ $^

build/longwave.o: longwave/src/longwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/disort_shortwave.o: shortwave/src/disort_shortwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/rayleigh.o: shortwave/src/rayleigh.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/shortwave.o: shortwave/src/shortwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/solar_flux.o: shortwave/src/solar_flux.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/RFM_voigt.o: gas-optics/src/RFM_voigt.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/cfcs.o: gas-optics/src/cfcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/collision_induced_absorption.o: gas-optics/src/collision_induced_absorption.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/cuda_kernels.o: gas-optics/src/cuda_kernels.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/gas_optics.o: gas-optics/src/gas_optics.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/kernel_utils.o: gas-optics/src/kernel_utils.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/kernels.o: gas-optics/src/kernels.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/launch.o: gas-optics/src/launch.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/molecules.o: gas-optics/src/molecules.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/ozone_continuum.o: gas-optics/src/ozone_continuum.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/parse_HITRAN_file.o: gas-optics/src/parse_HITRAN_file.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/spectral_bin.o: gas-optics/src/spectral_bin.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/tips2017.o: gas-optics/src/tips2017.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/water_vapor_continuum.o: gas-optics/src/water_vapor_continuum.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/argparse.o: utilities/src/argparse.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/curtis_godson.o: utilities/src/curtis_godson.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/device.o: utilities/src/device.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/optics.o: utilities/src/optics.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/parse_csv.o: utilities/src/parse_csv.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/spectral_grid.o: utilities/src/spectral_grid.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/utilities.o: utilities/src/utilities.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/verbosity.o: utilities/src/verbosity.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

# Cloud optics library.
build/clouds_lib.o: clouds/clouds_lib.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/hu_stamnes.o: clouds/hu_stamnes.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/ice_cloud_optics.o: clouds/ice_cloud_optics.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/incomplete_beta.o: clouds/incomplete_beta.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/netcdf_utils.o: clouds/netcdf_utils.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/optics_utils.o: clouds/optics_utils.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/stochastic_cloud.o: clouds/stochastic_clouds.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/libclouds.a: $(CLOUD_OBJECTS)
	ar rcs $@ $^

# Core library function tests.
build/test_harness.o: testing_harness/src/test_harness.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_longwave.o: longwave/test/test_longwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_disort_shortwave.o: shortwave/test/test_disort_shortwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_rayleigh.o: shortwave/test/test_rayleigh.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_shortwave.o: shortwave/test/test_shortwave.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_solar_flux.o: shortwave/test/test_solar_flux.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_RFM_voigt.o: gas-optics/test/test_RFM_voigt.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_cfcs.o: gas-optics/test/test_cfcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_collision_induced_absorption.o: gas-optics/test/test_collision_induced_absorption.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_cuda_kernels.o: gas-optics/test/test_cuda_kernels.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_gas_optics.o: gas-optics/test/test_gas_optics.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_kernel_utils.o: gas-optics/test/test_kernel_utils.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_kernels.o: gas-optics/test/test_kernels.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_launch.o: gas-optics/test/test_launch.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_molecules.o: gas-optics/test/test_molecules.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_ozone_continuum.o: gas-optics/test/test_ozone_continuum.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_parse_HITRAN_file.o: gas-optics/test/test_parse_HITRAN_file.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_spectral_bin.o: gas-optics/test/test_spectral_bin.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_tips2017.o: gas-optics/test/test_tips2017.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_water_vapor_continuum.o: gas-optics/test/test_water_vapor_continuum.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_argparse.o: utilities/test/test_argparse.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_curtis_godson.o: utilities/test/test_curtis_godson.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_device.o: utilities/test/test_device.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_optics.o: utilities/test/test_optics.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_parse_csv.o: utilities/test/test_parse_csv.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_spectral_grid.o: utilities/test/test_spectral_grid.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_utilities.o: utilities/test/test_utilities.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_verbosity.o: utilities/test/test_verbosity.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/test_longwave: build/test_longwave.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_disort_shortwave: build/test_disort_shortwave build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_rayleigh: build/test_rayleigh.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_shortwave: build/test_shortwave.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_solar_flux: build/test_solar_flux.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_RFM_voigt: build/test_RFM_voigt.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_cfcs: build/test_cfcs.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_collision_induced_absorption: build/test_collision_induced_absorption.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_cuda_kernels: build/test_cuda_kernels.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_gas_optics: build/test_gas_optics.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_kernel_utils: build/test_kernel_utils.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_kernels: build/test_kernels.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_launch: build/test_launch.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_molecules: build/test_molecules.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_ozone_continuum: build/test_ozone_continuum.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_parse_HITRAN_file: build/test_parse_HITRAN_file.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_spectral_bin: build/test_spectral_bin.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_tips2017: build/test_tips2017.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_water_vapor_continuum: build/test_water_vapor_continuum.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_argparse: build/test_argparse.o build/test_harness.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_curtis_godson: build/test_curtis_godson.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_device: build/test_device.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_optics: build/test_optics.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_parse_csv: build/test_parse_csv.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_spectral_grid: build/test_spectral_grid.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_utilities: build/test_utilities.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

build/test_verbosity: build/test_verbosity.o build/test_harness.o build/libgrtcode.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ -lm

# Applications.
# CIRC.
build/circ.o: circ/src/basic-circ-test.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/circ: framework/src/driver.c build/circ.o build/libclouds.a build/libgrtcode.a
	$(CC) $(CPPFLAGS) -Icirc/src $(CFLAGS) -o $@ $^ -lm

# Applications.
# ERA5.
build/era5.o: era5/src/era5.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/era5: framework/src/driver.c build/era5.o build/libclouds.a build/libgrtcode.a
	$(CC) $(CPPFLAGS) -Icirc/src $(CFLAGS) -o $@ $^ $(LDFLAGS) -lnetcdf -lm

# Applications.
# RFMIP-IRF
build/rfmip-irf.o: rfmip-irf/src/rfmip-irf.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

build/rfmip-irf: framework/src/driver.c build/rfmip-irf.o build/libclouds.a build/libgrtcode.a
	$(CC) $(CPPFLAGS) -Irfmip-irf/src $(CFLAGS) -o $@ $^ $(LDFLAGS) -lnetcdf -lm


clean:
	rm -f build/libgrtcode.a $(OBJECTS)
	rm -f build/libclouds.a $(CLOUD_OBJECTS)
	rm -f $(TEST_OBJECTS) $(TESTS)
	rm -f build/circ build/circ.o
	rm -f build/era5 build/era5.o
	rm -f build/rfmip-irf build/rfmip-irf.o
