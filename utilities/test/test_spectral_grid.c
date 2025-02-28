#include <stdlib.h>
#include "debug.h"
#include "device.h"
#include "floating_point_type.h"
#include "return_codes.h"
#include "spectral_grid.h"
#include "test_harness.h"
#include "utilities.h"


#define NUM_TESTS 13
#define DW 1.
#define N 3000
#define WN 3000.
#define W0 1.
#define TEST_GRID_SIZE 10


/*Helper function that initializes a spectral grid.*/
SpectralGrid_t initialize_spectral_grid()
{
    SpectralGrid_t grid = {
        .dw = DW,
        .n = N,
        .wn = WN,
        .w0 = W0
    };
    return grid;
}


/*Compare two identical grids.*/
int test_compare_spectral_grids_same()
{
    SpectralGrid_t one = initialize_spectral_grid();
    SpectralGrid_t two = initialize_spectral_grid();
    int result;
    rc_check(compare_spectral_grids(&one, &two, &result));
    if (!result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Compare grids with different lower bounds.*/
int test_compare_spectral_grids_different_start()
{
    SpectralGrid_t one = initialize_spectral_grid();
    SpectralGrid_t two = initialize_spectral_grid();
    two.w0 = one.w0 + 1.;
    int result;
    rc_check(compare_spectral_grids(&one, &two, &result));
    if (result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Compare grids with different upper bounds.*/
int test_compare_spectral_grids_different_stop()
{
    SpectralGrid_t one = initialize_spectral_grid();
    SpectralGrid_t two = initialize_spectral_grid();
    two.wn = one.wn - 1.;
    int result;
    rc_check(compare_spectral_grids(&one, &two, &result));
    if (result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Compare grids with different resolution.*/
int test_compare_spectral_grids_different_resolution()
{
    SpectralGrid_t one = initialize_spectral_grid();
    SpectralGrid_t two = initialize_spectral_grid();
    two.dw = 0.1*one.dw;
    int result;
    rc_check(compare_spectral_grids(&one, &two, &result));
    if (result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Create a grid.*/
int test_create_spectral_grid()
{
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, W0, WN, DW));
    SpectralGrid_t reference = initialize_spectral_grid();
    int result;
    rc_check(compare_spectral_grids(&grid, &reference, &result));
    if (!result)
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Create a grid with an invalid lower bound.*/
int test_create_spectral_grid_invalid_start()
{
    SpectralGrid_t grid;
    if (GRTCODE_RANGE_ERR != create_spectral_grid(&grid, -1., WN, DW))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Create a grid with an invalid upper bound.*/
int test_create_spectral_grid_invalid_stop()
{
    SpectralGrid_t grid;
    if (GRTCODE_RANGE_ERR != create_spectral_grid(&grid, W0, W0, DW))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Create a grid with an invalid resolution.*/
int test_create_spectral_grid_invalid_resolution()
{
    SpectralGrid_t grid;
    if (GRTCODE_RANGE_ERR != create_spectral_grid(&grid, W0, W0, -1.))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Find the index of a point on the grid.*/
int test_grid_point_index()
{
    SpectralGrid_t grid = initialize_spectral_grid();
    double w = WN - W0;
    uint64_t index;
    rc_check(grid_point_index(grid, w, &index));
    if (index != (uint64_t)(((int)(1./DW))*(w - W0)))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Look for a grid point that is outside of the range covered by the grid.*/
int test_grid_point_index_outside_of_grid()
{
    SpectralGrid_t grid = initialize_spectral_grid();
    double w = WN + 50;
    uint64_t index;
    if (GRTCODE_RANGE_ERR != grid_point_index(grid, w, &index))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Look for a grid point that is not on the grid.*/
int test_grid_point_index_not_on_grid()
{
    SpectralGrid_t grid = initialize_spectral_grid();
    double w = W0 + 0.5*DW;
    uint64_t index;
    if (GRTCODE_VALUE_ERR != grid_point_index(grid, w, &index))
    {
        return GRTCODE_VALUE_ERR;
    }
    return GRTCODE_SUCCESS;
}


/*Create an array of grid points.*/
int test_grid_points()
{
    double w0 = 1.;
    double wn = 10000.;
    double dw = 0.1;
    SpectralGrid_t grid;
    rc_check(create_spectral_grid(&grid, w0, wn, dw));
    Device_t device;
    rc_check(create_device(&device, NULL));
    fp_t * buffer;
    rc_check(grid_points(grid, &buffer, device));
    int returncode;
    if (buffer[0] == w0 && buffer[grid.n - 1] == wn)
    {
        returncode = GRTCODE_SUCCESS;
    }
    else
    {
        returncode = GRTCODE_VALUE_ERR;
    }
    gfree(buffer, device);
    return returncode;
}


/*Interpolate data on to a grid.*/
int test_interpolate_to_grid()
{
    SpectralGrid_t grid = initialize_spectral_grid();
    fp_t x[TEST_GRID_SIZE];
    fp_t y[TEST_GRID_SIZE];
    int i;
    for (i=0; i<TEST_GRID_SIZE; ++i)
    {
        x[i] = W0 + (fp_t)i + 0.25;
        y[i] = (fp_t)i;
    }
    fp_t * newy = (fp_t *)malloc(sizeof(*newy)*grid.n);
    Sample1d_t interp = linear_sample;
    Sample1d_t extrap = constant_extrapolation;
    rc_check(interpolate_to_grid(grid, x, y, TEST_GRID_SIZE, newy, interp, extrap));
    free(newy);
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_compare_spectral_grids_same,
            "test_compare_spectral_grids_same",
            GRTCODE_SUCCESS
        },
        {
            test_compare_spectral_grids_different_start,
            "test_compare_spectral_grids_different_start",
            GRTCODE_SUCCESS
        },
        {
            test_compare_spectral_grids_different_stop,
            "test_compare_spectral_grids_different_stop",
            GRTCODE_SUCCESS
        },
        {
            test_compare_spectral_grids_different_resolution,
            "test_compare_spectral_grids_different_resolution",
            GRTCODE_SUCCESS
        },
        {
            test_create_spectral_grid,
            "test_create_spectral_grid",
            GRTCODE_SUCCESS
        },
        {
            test_create_spectral_grid_invalid_start,
            "test_create_spectral_grid_invalid_start",
            GRTCODE_SUCCESS
        },
        {
            test_create_spectral_grid_invalid_stop,
            "test_create_spectral_grid_invalid_stop",
            GRTCODE_SUCCESS
        },
        {
            test_create_spectral_grid_invalid_resolution,
            "test_create_spectral_grid_invalid_resolution",
            GRTCODE_SUCCESS
        },
        {
            test_grid_point_index,
            "test_grid_point_index",
            GRTCODE_SUCCESS
        },
        {
            test_grid_point_index_outside_of_grid,
            "test_grid_point_index_outside_of_grid",
            GRTCODE_SUCCESS
        },
        {
            test_grid_point_index_not_on_grid,
            "test_grid_point_index_not_on_grid",
            GRTCODE_SUCCESS
        },
        {
            test_grid_points,
            "test_grid_points",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate_to_grid,
            "test_interpolate_to_grid",
            GRTCODE_SUCCESS
        },
    };
    return test_harness("test_spectral_grid", NUM_TESTS, tests);
}
