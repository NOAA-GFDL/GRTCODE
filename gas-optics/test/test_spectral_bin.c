#include "grtcode_utilities.h"
#include "spectral_bin.h"
#include "spectral_bin-internal.h"
#include "test_harness.h"


#define NUM_TESTS 2


int test_create_spectral_bins()
{
    SpectralBins_t bins;
    int num_layers = 5;
    double w0 = 1.;
    uint64_t n = 100;
    double wres = 0.1;
    double bin_width = 1.;
    Device_t device;
    rc_check(create_spectral_bins(&bins, num_layers, w0, n, wres, bin_width, device));
    return GRTCODE_SUCCESS;
}


int test_destroy_spectral_bins()
{
    SpectralBins_t bins;
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_create_spectral_bins,
            "test_create_spectral_bins",
            GRTCODE_SUCCESS
        },
        {
            test_destroy_spectral_bins,
            "test_destroy_spectral_bins",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_spectral_bins", NUM_TESTS, tests);
}
