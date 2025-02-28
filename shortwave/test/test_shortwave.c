#include <stdint.h>
#include <stdlib.h>
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "shortwave.h"
#include "test_harness.h"


#define NUM_TESTS 4


typedef struct Atmosphere
{
    Device_t device;
    int num_levels;
    fp_t mu_dir;
    fp_t mu_dif;
    fp_t *surface_albedo_dir;
    fp_t *surface_albedo_dif;
    fp_t total_solar_irradiance;
    fp_t *solar_flux;
    fp_t *flux_up;
    fp_t *flux_down;
} Atmosphere_t;


static int create_atmosphere(Atmosphere_t * const atmos, SpectralGrid_t const * const grid,
                             int const num_levels, Device_t const device)
{
    atmos->device = device;
    atmos->num_levels = num_levels;
    atmos->surface_albedo_dir = NULL;
    gmalloc(atmos->surface_albedo_dir, grid->n, device);
    atmos->surface_albedo_dif = NULL;
    gmalloc(atmos->surface_albedo_dif, grid->n, device);
    atmos->solar_flux = NULL;
    gmalloc(atmos->solar_flux, grid->n, device);
    atmos->flux_up = NULL;
    gmalloc(atmos->flux_up, grid->n*num_levels, device);
    atmos->flux_down = NULL;
    gmalloc(atmos->flux_down, grid->n*num_levels, device);
    return GRTCODE_SUCCESS;
}


static int destroy_atmosphere(Atmosphere_t * const atmos)
{
    gfree(atmos->surface_albedo_dir, atmos->device);
    gfree(atmos->surface_albedo_dif, atmos->device);
    gfree(atmos->solar_flux, atmos->device);
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
            optics->g[i*grid->n+j] = 0.85*(rand()/((double)RAND_MAX));
            optics->omega[i*grid->n+j] = 0.85*(rand()/((double)RAND_MAX));
            optics->tau[i*grid->n+j] = rand()/((double)RAND_MAX);
        }
    }
    catch(create_atmosphere(atmos, grid, num_levels, device));
    for (j=0; j<grid->n; ++j)
    {
        atmos->surface_albedo_dir[j] = 0.35;
        atmos->surface_albedo_dif[j] = 0.35;
        atmos->solar_flux[j] = 0.99*rand()/((double)RAND_MAX);
    }
    atmos->mu_dir = 0.8;
    atmos->mu_dif = 0.5;
    atmos->total_solar_irradiance = 1350.;
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
    Shortwave_t sw;
    catch(create_shortwave(&sw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_sw_fluxes(&sw, &optics, atmos.mu_dir, atmos.mu_dif,
                              atmos.surface_albedo_dir, atmos.surface_albedo_dif,
                              atmos.total_solar_irradiance, atmos.solar_flux, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_shortwave(&sw));
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
    Shortwave_t sw;
    catch(create_shortwave(&sw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_sw_fluxes(&sw, &optics, atmos.mu_dir, atmos.mu_dif,
                              atmos.surface_albedo_dir, atmos.surface_albedo_dif,
                              atmos.total_solar_irradiance, atmos.solar_flux, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_shortwave(&sw));
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
    Shortwave_t sw;
    catch(create_shortwave(&sw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_sw_fluxes(&sw, &optics, atmos.mu_dir, atmos.mu_dif,
                              atmos.surface_albedo_dir, atmos.surface_albedo_dif,
                              atmos.total_solar_irradiance, atmos.solar_flux, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_shortwave(&sw));
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
    Shortwave_t sw;
    catch(create_shortwave(&sw, atmos.num_levels, &grid, &atmos.device));
    catch(calculate_sw_fluxes(&sw, &optics, atmos.mu_dir, atmos.mu_dif,
                              atmos.surface_albedo_dir, atmos.surface_albedo_dif,
                              atmos.total_solar_irradiance, atmos.solar_flux, atmos.flux_up,
                              atmos.flux_down));
    catch(destroy_shortwave(&sw));
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
    return test_harness("test_shortwave", NUM_TESTS, tests);
}
