#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "return_codes.h"
#include "test_harness.h"
#include "utilities.h"


#define NUM_TESTS 17


int test_activate()
{
    uint64_t bit_field = 0;
    int index = 5;
    rc_check(activate(&bit_field, index));
    uint64_t result = pow(2, index - 1);
    if (bit_field != result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_angstrom_exponent()
{
    fp_t tau1 = 3.;
    fp_t tau2 = 5.5;
    fp_t lambda1 = 7.;
    fp_t lambda2 = 9.1;
    fp_t result = -2.31028339473;
    if (result != angstrom_exponent(tau1, tau2, lambda1, lambda2))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_angstrom_exponent_sample()
{
    fp_t x[2] = {5., 25.};
    fp_t y[2] = {0.4567, 2.9801};
    fp_t newx[32];
    fp_t newy[32];
    size_t n = 32;
    rc_check(angstrom_exponent_sample(x, y, newx, newy, n));
    return GRTCODE_SUCCESS;
}


int test_constant_extrapolation()
{
    fp_t x[32];
    fp_t y[32];
    fp_t newx[32];
    fp_t newy[32];
    size_t n = 32;
    rc_check(constant_extrapolation(x, y, newx, newy, n));
    return GRTCODE_SUCCESS;
}


int test_copy_str()
{
    char destination[32];
    char source[32];
    size_t length = 32;
    snprintf(source, length, "This is a test.");
    rc_check(copy_str(destination, source, length));
    return GRTCODE_SUCCESS;
}


int test_free_ptr()
{
    void * p = malloc(128);
    rc_check(free_ptr(&p));
    return GRTCODE_SUCCESS;
}


int test_integrate2()
{
    fp_t x[32];
    fp_t y[32];
    size_t n = 32;
    fp_t s;
    Area1d_t area = trapezoid;
    rc_check(integrate2(x, y, n, &s, area));
    return GRTCODE_SUCCESS;
}


int test_interpolate2()
{
    fp_t x[32];
    fp_t y[32];
    size_t n = 32;
    fp_t newx[32];
    fp_t newy[32];
    size_t newn = 32;
    Sample1d_t interp = linear_sample;
    Sample1d_t extrap = constant_extrapolation;
    rc_check(interpolate2(x, y, n, newx, newy, newn, interp, extrap));
    return GRTCODE_SUCCESS;
}


int test_is_active()
{
    uint64_t bit_field = 16;
    int index = 5;
    if (!is_active(bit_field, index))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_linear_sample()
{
    fp_t x[32];
    fp_t y[32];
    fp_t newx[32];
    fp_t newy[32];
    size_t n = 32;
    rc_check(linear_sample(x, y, newx, newy, n));
    return GRTCODE_SUCCESS;
}


int test_malloc_ptr()
{
    void * p;
    rc_check(malloc_ptr(&p, 128));
    return GRTCODE_SUCCESS;
}


int test_monotonically_increasing()
{
    fp_t x[32];
    size_t n = 32;
    if (!monotonically_increasing(x, n))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_open_file()
{
    FILE * file;
    char * path = "foo";
    char * mode = "r";
    rc_check(open_file(&file, path, mode));
    return GRTCODE_SUCCESS;
}


int test_to_double()
{
    char * s = "1.";
    double d;
    rc_check(to_double(s, &d));
    if (d != 1.)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_to_fp_t()
{
    double d = 1.;
    fp_t f;
    rc_check(to_fp_t(d, &f));
    if (f != (fp_t)1.)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_to_int()
{
    char * s = "55";
    int i;
    rc_check(to_int(s, &i));
    if (i != 55)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_trapezoid()
{
    fp_t x[32];
    fp_t y[32];
    fp_t result = trapezoid(x, y);
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_activate,
            "test_activate",
            GRTCODE_SUCCESS
        },
        {
            test_angstrom_exponent,
            "test_angstrom_exponent",
            GRTCODE_SUCCESS
        },
        {
            test_angstrom_exponent_sample,
            "test_angstrom_exponent_sample",
            GRTCODE_SUCCESS
        },
        {
            test_constant_extrapolation,
            "test_constant_extrapolation",
            GRTCODE_SUCCESS
        },
        {
            test_copy_str,
            "test_copy_str",
            GRTCODE_SUCCESS
        },
        {
            test_free_ptr,
            "test_free_ptr",
            GRTCODE_SUCCESS
        },
        {
            test_integrate2,
            "test_integrate2",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate2,
            "test_interpolate2",
            GRTCODE_SUCCESS
        },
        {
            test_is_active,
            "test_is_active",
            GRTCODE_SUCCESS
        },
        {
            test_linear_sample,
            "test_linear_sample",
            GRTCODE_SUCCESS
        },
        {
            test_malloc_ptr,
            "test_malloc_ptr",
            GRTCODE_SUCCESS
        },
        {
            test_monotonically_increasing,
            "test_monotonically_increasing",
            GRTCODE_SUCCESS
        },
        {
            test_open_file,
            "test_open_file",
            GRTCODE_SUCCESS
        },
        {
            test_to_double,
            "test_to_double",
            GRTCODE_SUCCESS
        },
        {
            test_to_fp_t,
            "test_to_fp_t",
            GRTCODE_SUCCESS
        },
        {
            test_to_int,
            "test_to_int",
            GRTCODE_SUCCESS
        },
        {
            test_trapezoid,
            "test_trapezoid",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_utilities", NUM_TESTS, tests);
}
