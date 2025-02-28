#include <stdint.h>
#include "grtcode_utilities.h"
#include "kernel_utils.h"
#include "test_harness.h"


#define NUM_TESTS 3


int test_bracket()
{
    uint64_t array_size = 10;
    fp_t array[10] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
    fp_t value = 5.5;
    uint64_t left;
    uint64_t right;
    rc_check(bracket(array_size, array, value, &left, &right));
    if (left != 4 || right != 5)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_bin_quad_interp()
{
    fp_t x[10];
    fp_t y[10];
    uint64_t left = 0;
    uint64_t right = 9;
    fp_t w0 = 1.;
    fp_t dw = 0.1;
    fp_t tau[10];
    rc_check(bin_quad_interp(x, y, left, right, w0, dw, tau));
    return GRTCODE_SUCCESS;
}


int test_bin_no_interp()
{
    uint64_t left = 0;
    uint64_t right = 9;
    fp_t taub[10];
    fp_t tau[10];
    rc_check(bin_no_interp(left, right, taub, tau));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_bracket,
            "test_bracket",
            GRTCODE_SUCCESS
        },
        {
            test_bin_quad_interp,
            "test_bin_quad_interp",
            GRTCODE_SUCCESS
        },
        {
            test_bin_no_interp,
            "test_bin_no_interp",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_kernel_utils", NUM_TESTS, tests);
}
