#include "grtcode_utilities.h"
#include "ozone_continuum.h"
#include "ozone_continuum-internal.h"
#include "test_harness.h"


#define NUM_TESTS 2


int test_get_ozone_continuum_coefs()
{
    OzoneContinuumCoefs_t cc;
    char * o3_ctm_file = "";
    SpectralGrid_t grid;
    Device_t device;
    rc_check(get_ozone_continuum_coefs(&cc, o3_ctm_file, grid, device));
    return GRTCODE_SUCCESS;
}


int test_free_ozone_continuum_coefs()
{
    OzoneContinuumCoefs_t cc;
    rc_check(free_ozone_continuum_coefs(&cc));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_get_ozone_continuum_coefs,
            "test_get_ozone_continuum_coefs",
            GRTCODE_SUCCESS
        },
        {
            test_free_ozone_continuum_coefs,
            "test_free_ozone_continuum_coefs",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_ozone_continuum", NUM_TESTS, tests);
}
