#include "grtcode_utilities.h"
#include "test_harness.h"
#include "water_vapor_continuum.h"


#define NUM_TESTS 2


int test_get_water_vapor_continuum_coefs()
{
    WaterVaporContinuumCoefs_t cc;
    char * h2o_ctm_dir = "";
    SpectralGrid_t grid;
    Device_t device;
    rc_check(get_water_vapor_continuum_coefs(&cc, h2o_ctm_dir, grid, device));
    return GRTCODE_SUCCESS;
}


int test_free_water_vapor_continuum_coefs()
{
    WaterVaporContinuumCoefs_t cc;
    rc_check(free_water_vapor_continuum_coefs(&cc));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_get_water_vapor_continuum_coefs,
            "test_get_water_vapor_continuum_coefs",
            GRTCODE_SUCCESS
        },
        {
            test_free_water_vapor_continuum_coefs,
            "test_free_water_vapor_continuum_coefs",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_water_vapor_continuum", NUM_TESTS, tests);
}
