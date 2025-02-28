#include <stdint.h>

#include "floating_point_type.h"
#include "optics.h"
#include "return_codes.h"
#include "test_harness.h"


#define NUM_TESTS 6


/*Add optical properties together.*/
int test_add_optics()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t one;
    int num_layers = 5;
    rc_check(create_optics(&one, num_layers, &grid, &device));
    Optics_t two;
    rc_check(create_optics(&two, num_layers, &grid, &device));
    Optics_t const * optics_array[2] = {&one, &two};
    int num_optics = 2;
    Optics_t result;
    rc_check(add_optics(optics_array, num_optics, &result));
    rc_check(destroy_optics(&one));
    rc_check(destroy_optics(&two));
    return GRTCODE_SUCCESS;
}


int test_create_optics()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t gas;
    int num_layers = 5;
    rc_check(create_optics(&gas, num_layers, &grid, &device));
    int returncode_value;
    if (gas.num_layers != num_layers || gas.grid.n != grid.n || gas.device != device)
    {
        returncode_value = GRTCODE_VALUE_ERR;
    }
    else
    {
        returncode_value = GRTCODE_SUCCESS;
    }
    rc_check(destroy_optics(&gas));
    return returncode_value;
}


int test_destroy_optics()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t optics;
    int num_layers = 5;
    rc_check(create_optics(&optics, num_layers, &grid, &device));
    rc_check(destroy_optics(&optics));
    return GRTCODE_SUCCESS;
}


int test_optics_compatible()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t one;
    int num_layers = 5;
    rc_check(create_optics(&one, num_layers, &grid, &device));
    Optics_t two;
    rc_check(create_optics(&two, num_layers, &grid, &device));
    int result;
    rc_check(optics_compatible(&one, &two, &result));
    int returncode_value;
    if (!result)
    {
        returncode_value = GRTCODE_VALUE_ERR;
    }
    else
    {
        returncode_value = GRTCODE_SUCCESS;
    }
    rc_check(destroy_optics(&one));
    rc_check(destroy_optics(&two));
    return returncode_value;
}


int test_sample_optics()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t destination;
    int num_layers = 5;
    rc_check(create_optics(&destination, num_layers, &grid, &device));
    Optics_t source;
    rc_check(create_optics(&source, num_layers, &grid, &device));
    rc_check(sample_optics(&destination, &source, &w0, &wn));
    rc_check(destroy_optics(&destination));
    rc_check(destroy_optics(&source));
    return GRTCODE_SUCCESS;
}


int test_update_optics()
{
    Device_t device;
    rc_check(create_device(&device, NULL));
    SpectralGrid_t grid;
    double w0 = 1.;
    double wn = 100.;
    double dw = 0.1;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Optics_t optics;
    int num_layers = 5;
    rc_check(create_optics(&optics, num_layers, &grid, &device));
    fp_t * tau = (fp_t *)malloc(sizeof(*tau)*grid.n);
    fp_t * omega = (fp_t *)malloc(sizeof(*omega)*grid.n);
    fp_t * g = (fp_t *)malloc(sizeof(*g)*grid.n);
    rc_check(update_optics(&optics, tau, omega, g));
    rc_check(destroy_optics(&optics));
    free(tau);
    free(omega);
    free(g);
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_add_optics,
            "test_add_optics",
            GRTCODE_SUCCESS
        },
        {
            test_create_optics,
            "test_create_optics",
            GRTCODE_SUCCESS
        },
        {
            test_destroy_optics,
            "test_destroy_optics",
            GRTCODE_SUCCESS
        },
        {
            test_optics_compatible,
            "test_optics_compatible",
            GRTCODE_SUCCESS
        },
        {
            test_sample_optics,
            "test_sample_optics",
            GRTCODE_SUCCESS
        },
        {
            test_update_optics,
            "test_update_optics",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_optics", NUM_TESTS, tests);
}
