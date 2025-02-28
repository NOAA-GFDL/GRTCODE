#include "cfcs.h"
#include "cfcs-internal.h"
#include "grtcode_utilities.h"
#include "test_harness.h"


#define NUM_TESTS 2


int test_get_cfc_cross_sections()
{
    Device_t device;
    SpectralGrid_t grid;
    int id = 0;
    char * filepath = "";
    CfcCrossSection_t cross_section;
    rc_check(get_cfc_cross_sections(&cross_section, id, filepath, grid, device));
    return GRTCODE_SUCCESS;
}


int test_free_cfc_cross_sections()
{
    CfcCrossSection_t cross_section;
    rc_check(free_cfc_cross_sections(&cross_section));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_get_cfc_cross_sections,
            "test_get_cfc_cross_sections",
            GRTCODE_SUCCESS
        },
        {
            test_free_cfc_cross_sections,
            "test_free_cfc_cross_sections",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_cfcs", NUM_TESTS, tests);
}
