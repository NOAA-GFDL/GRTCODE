#ifndef TEST_HARNESS_H_
#define TEST_HARNESS_H_


#include <stdlib.h>
#include "device.h"
#include "floating_point_type.h"


#define rc_check(function) { \
    int rc = function; \
    if (rc != GRTCODE_SUCCESS) \
    {  \
        return GRTCODE_VALUE_ERR; \
    } \
}


static Device_t device = HOST_ONLY;


/*Check values in an array.*/
int check_array(
    fp_t const * const actual, /*Array of actual values.*/
    fp_t const * const expected, /*Array of expected values.*/
    size_t const size, /*Size of the input arrays.*/
    fp_t const relative_tolerance, /*Relative Tolerance.*/
    fp_t const absolute_tolerance /*Use absolute differences or percent differences?*/
);


/*Compare two floating point values.*/
int check_floating_point(
    fp_t const actual, /*Actual value.*/
    fp_t const expected, /*Expected value.*/
    fp_t const relative_tolerance, /*Relative Tolerance.*/
    fp_t const absolute_tolerance /*Use absolute differences or percent differences?*/
);


/*Helper structure so tests can be looped through.*/
typedef struct Test
{
    int (* test_function_pointer)(void);
    char * name;
    int returncode;
} Test_t;


/*Testing harness.*/
int test_harness(
    char const * const title, /*Name of the test suite.*/
    int const num_tests, /*Size of the input tests array.*/
    Test_t const * const tests /*Array of tests.*/
);


#endif
