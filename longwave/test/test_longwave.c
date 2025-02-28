#include <stdint.h>
#include <stdlib.h>
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "longwave.h"
#include "test_harness.h"


#define NUM_TESTS 4


typedef struct Atmosphere
{
    Device_t device;
    int num_levels;
    fp_t *layer_temperature;
    fp_t *level_temperature;
    fp_t surface_temperature;
    fp_t *emissivity;
    fp_t *flux_up;
    fp_t *flux_down;
} Atmosphere_t;


static int create_atmosphere(Atmosphere_t * const atmos, SpectralGrid_t const * const grid,
                             int const num_levels, Device_t const device)
{
    int num_layers = num_levels - 1;
    atmos->device = device;
    atmos->num_levels = num_levels;
    atmos->layer_temperature = NULL;
    gmalloc(atmos->layer_temperature, num_layers, device);
    atmos->level_temperature = NULL;
    gmalloc(atmos->level_temperature, num_levels, device);
    atmos->emissivity = NULL;
    gmalloc(atmos->emissivity, grid->n, device);
    atmos->flux_up = NULL;
    gmalloc(atmos->flux_up, grid->n*num_levels, device);
    atmos->flux_down = NULL;
    gmalloc(atmos->flux_down, grid->n*num_levels, device);
    return GRTCODE_SUCCESS;
}


static int destroy_atmosphere(Atmosphere_t * const atmos)
{
    gfree(atmos->layer_temperature, atmos->device);
    gfree(atmos->level_temperature, atmos->device);
    gfree(atmos->emissivity, atmos->device);
    gfree(atmos->flux_up, atmos->device);
    gfree(atmos->flux_down, atmos->device);
    return GRTCODE_SUCCESS;
}


static int setup(SpectralGrid_t * const grid, Optics_t * const optics,
                 Atmosphere_t * const atmos)
{
    Device_t const device = HOST_ONLY;
    double const w0 = MIN_WAVENUMBER;
    double const wn = MAX_WAVENUMBER;
    double const dw = MAX_RESOLUTION;
    int const num_levels = MIN_NUM_LEVELS;
    int const num_layers = num_levels - 1;
    catch(create_spectral_grid(grid, w0, wn, dw));
    catch(create_optics(optics, num_layers, grid, &device));
    int i;
    uint64_t j;
    for (i=0; i<num_layers; ++i)
    {
        for (j=0; j<grid->n; ++j)
        {
            optics->tau[i*grid->n+j] = rand()/((double)RAND_MAX);
        }
    }
    catch(create_atmosphere(atmos, grid, num_levels, device));
    for (i=0; i<num_layers; ++i)
    {
        atmos->layer_temperature[i] = 250. + 1.*i;
        atmos->level_temperature[i] = atmos->layer_temperature[i] - 0.5;
    }
    atmos->level_temperature[num_layers] = atmos->layer_temperature[num_layers-1] + 0.5;
    atmos->surface_temperature = atmos->level_temperature[num_layers] + 2.;
    for (j=0; j<grid->n; ++j)
    {
        atmos->emissivity[j] = 0.98;
    }
    return GRTCODE_SUCCESS;
}


static int breakdown(Optics_t * const optics, Atmosphere_t * const atmos)
{
    catch(destroy_atmosphere(atmos));
    catch(destroy_optics(optics));
    return GRTCODE_SUCCESS;
}


/*Do a simple run.*/
int test_simple()
{
    SpectralGrid_t grid;
    Optics_t optics;
    Atmosphere_t atmos;
    catch(setup(&grid, &optics, &atmos));
    Longwave_t lw;
    catch(create_longwave(&lw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_lw_fluxes(&lw, &optics, atmos.surface_temperature, atmos.layer_temperature,
                              atmos.level_temperature, atmos.emissivity, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_longwave(&lw));
    catch(breakdown(&optics, &atmos));
    return GRTCODE_SUCCESS;
}


/*Very optically thick layers.*/
int test_optically_thick()
{
    SpectralGrid_t grid;
    Optics_t optics;
    Atmosphere_t atmos;
    catch(setup(&grid, &optics, &atmos));
    int i;
    uint64_t j;
    for (i=0; i<(atmos.num_levels-1); ++i)
    {
        for (j=0; j<grid.n; ++j)
        {
            /*Big enough to overflow exp(l.tau).*/
            optics.tau[i*grid.n+j] = 1.e12;
        }
    }
    Longwave_t lw;
    catch(create_longwave(&lw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_lw_fluxes(&lw, &optics, atmos.surface_temperature, atmos.layer_temperature,
                              atmos.level_temperature, atmos.emissivity, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_longwave(&lw));
    catch(breakdown(&optics, &atmos));
    return GRTCODE_SUCCESS;
}


/*Very optically thin layers.*/
int test_optically_thin()
{
    SpectralGrid_t grid;
    Optics_t optics;
    Atmosphere_t atmos;
    catch(setup(&grid, &optics, &atmos));
    int i;
    uint64_t j;
    for (i=0; i<(atmos.num_levels-1); ++i)
    {
        for (j=0; j<grid.n; ++j)
        {
            /*Small enough to underflow exp(l.tau).*/
            optics.tau[i*grid.n+j] = 1.e-12;
        }
    }
    Longwave_t lw;
    catch(create_longwave(&lw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_lw_fluxes(&lw, &optics, atmos.surface_temperature, atmos.layer_temperature,
                              atmos.level_temperature, atmos.emissivity, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_longwave(&lw));
    catch(breakdown(&optics, &atmos));
    return GRTCODE_SUCCESS;
}


/*Strong absorption.*/
int test_strong_absorption()
{
    SpectralGrid_t grid;
    Optics_t optics;
    Atmosphere_t atmos;
    catch(setup(&grid, &optics, &atmos));
    int nbytes = sizeof(*optics.tau);
    fp_t val = 0.;
    if (nbytes == sizeof(float))
    {
        val = log(FLT_MAX) - 12.;
    }
    else if (nbytes == sizeof(double))
    {
        val = log(DBL_MAX) - 12.;
    }
    int i;
    uint64_t j;
    for (i=0; i<(atmos.num_levels-1); ++i)
    {
        for (j=0; j<grid.n; ++j)
        {
            optics.tau[i*grid.n+j] = val;
        }
    }
    Longwave_t lw;
    catch(create_longwave(&lw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_lw_fluxes(&lw, &optics, atmos.surface_temperature, atmos.layer_temperature,
                              atmos.level_temperature, atmos.emissivity, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_longwave(&lw));
    catch(breakdown(&optics, &atmos));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_simple,
            "test_simple",
            GRTCODE_SUCCESS
        },
        {
            test_optically_thick,
            "test_optically_thick",
            GRTCODE_SUCCESS
        },
        {
            test_optically_thin,
            "test_optically_thin",
            GRTCODE_SUCCESS
        },
        {
            test_strong_absorption,
            "test_strong_absorption",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_longwave", NUM_TESTS, tests);
}
