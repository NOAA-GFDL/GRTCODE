#include "grtcode_utilities.h"
#include "molecules.h"
#include "parse_HITRAN_file.h"
#include "parse_HITRAN_file-internal.h"
#include "test_harness.h"


#define NUM_TESTS 2


/*Free memory used for molecular line parameters.*/
int test_free_line_parameters()
{
    LineParams_t line_params;
    rc_check(free_line_params(&line_params));
    return GRTCODE_SUCCESS;
}


int test_parse_hitran_file()
{
    LineParams_t line_params;
    char * filepath = "";
    int mol_id = H2O;
    double w0 = 1.;
    double wn = 1000.;
    Device_t device;
    rc_check(parse_hitran_file(&line_params, filepath, mol_id, w0, wn, device));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_free_line_parameters,
            "test_free_line_parameters",
            GRTCODE_SUCCESS
        },
        {
            test_parse_hitran_file,
            "test_parse_hitran_file",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_parse_HITRAN_file", NUM_TESTS, tests);
}
