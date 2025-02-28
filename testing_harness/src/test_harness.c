#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "floating_point_type.h"
#include "test_harness.h"


/*Check values in an array.*/
int check_array(fp_t const * const actual, fp_t const * const expected, size_t const size,
                fp_t const relative_tolerance, fp_t const absolute_tolerance)
{
    int failures = 0;
    size_t i;
    for (i=0; i<size; ++i)
    {
        if (!check_floating_point(actual[i], expected[i], relative_tolerance, absolute_tolerance))
        {
            fprintf(stderr, "\tLarger than expect difference found at index [%d]\n", (int)i);
            failures++;
        }
    }
    return failures;
}


/*Return the maximum of two floating point values.*/
fp_t maximum(fp_t a, fp_t b)
{
    return a > b ? a : b;
}


/*Check if two floating point numbers are approximately equal.*/
int check_floating_point(fp_t const actual, fp_t const expected,
                         fp_t const relative_tolerance, fp_t const absolute_tolerance)
{
    fp_t value = maximum(fabs(actual), fabs(expected));
    fp_t tolerance = maximum(relative_tolerance*value, absolute_tolerance);
    if (fabs(actual - expected) <= tolerance)
    {
        return 1;
    }
    else
    {
        fprintf(stderr, "Error: actual: %e,  expected: %e,  difference: %e \n",
                actual, expected, actual - expected);
        return 0;
    }
}


/*Testing harness.*/
int test_harness(char const * const title, int const num_tests, Test_t const * const tests)
{
    printf("Running %s tests.\n", title);
    int passed = 0;
    int i;
    for (i=0; i<num_tests; ++i)
    {
        printf("Running test %s:", tests[i].name);
        if (tests[i].test_function_pointer() == tests[i].returncode)
        {
            printf(" passed.\n");
            passed++;
        }
        else
        {
            printf(" failed.\n");
        }
    }
    printf("%d/%d tests passed.\n", passed, num_tests);
    if (passed == num_tests)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
