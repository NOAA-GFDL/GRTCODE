#include <string.h>
#include "return_codes.h"
#include "test_harness.h"
#include "verbosity.h"
#include "verbosity-internal.h"


#define NUM_TESTS 6


int test_append_to_error_buffer()
{
    append_to_error_buffer("test string.");
    return GRTCODE_SUCCESS;
}


int test_copy_error_buffer()
{
    char buffer[4096];
    copy_error_buffer(buffer, 4096);
    return GRTCODE_SUCCESS;
}


int test_grtcode_errstr()
{
    char buffer[4096];
    rc_check(grtcode_errstr(GRTCODE_SUCCESS, buffer, 4096));
    return GRTCODE_SUCCESS;
}


int test_grtcode_set_verbosity()
{
    grtcode_set_verbosity(1);
    return GRTCODE_SUCCESS;
}


int test_grtcode_verbosity()
{
    int level = 2;
    grtcode_set_verbosity(level);
    int returned_level = grtcode_verbosity();
    if (level != returned_level)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


int test_reset_error_buffer()
{
    reset_error_buffer();
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_append_to_error_buffer,
            "test_append_to_error_buffer",
            GRTCODE_SUCCESS
        },
        {
            test_copy_error_buffer,
            "test_copy_error_buffer",
            GRTCODE_SUCCESS
        },
        {
            test_grtcode_errstr,
            "test_grtcode_errstr",
            GRTCODE_SUCCESS
        },
        {
            test_grtcode_set_verbosity,
            "test_grtcode_set_verbosity",
            GRTCODE_SUCCESS
        },
        {
            test_grtcode_verbosity,
            "test_grtcode_verbosity",
            GRTCODE_SUCCESS
        },
        {
            test_reset_error_buffer,
            "test_reset_error_buffer",
            GRTCODE_SUCCESS
        },
    };
    return test_harness("test_verbosity", NUM_TESTS, tests);
}
