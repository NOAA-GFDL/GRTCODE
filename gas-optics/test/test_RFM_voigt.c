#include "line_shape.h"
#include "return_codes.h"
#include "RFM_voigt.h"
#include "test_harness.h"


#define NUM_TESTS 1


int test_rfm_voigt_line_shape()
{
    LineShapeInputs_t inputs;
    int n = 1000;
    fp_t line_shape_function[1000];
    rc_check(rfm_voigt_line_shape(inputs, line_shape_function));
    fp_t reference[1000];
    fp_t tolerance = 1.e-6;
    rc_check(check_array(reference, line_shape_function, n, tolerance, 1));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_rfm_voigt_line_shape,
            "test_rfm_voigt_line_shape",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_RFM_voigt", NUM_TESTS, tests);
}
