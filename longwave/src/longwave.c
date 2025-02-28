/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include <stdint.h>
#include <string.h>
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "longwave.h"


/*Reserve memory for the longwave.*/
EXTERN int create_longwave(Longwave_t * const lw, int const num_levels,
                           SpectralGrid_t const * const grid, Device_t const * const device)
{
    not_null(lw);
    not_null(grid);
    not_null(device);
    in_range(num_levels, MIN_NUM_LEVELS, MAX_NUM_LEVELS);
    lw->num_levels = num_levels;
    lw->grid = *grid;
    lw->device = *device;
    if (lw->device != HOST_ONLY)
    {
        gmalloc(lw->layer_temperature, num_levels-1, lw->device);
        gmalloc(lw->level_temperature, num_levels, lw->device);
        gmalloc(lw->emissivity, lw->grid.n, lw->device);
        gmalloc(lw->flux_up, lw->grid.n*num_levels, lw->device);
        gmalloc(lw->flux_down, lw->grid.n*num_levels, lw->device);
    }
    return GRTCODE_SUCCESS;
}


/*Free memory for the longwave.*/
EXTERN int destroy_longwave(Longwave_t * const lw)
{
    not_null(lw);
    if (lw->device != HOST_ONLY)
    {
        gfree(lw->layer_temperature, lw->device);
        gfree(lw->level_temperature, lw->device);
        gfree(lw->emissivity, lw->device);
        gfree(lw->flux_up, lw->device);
        gfree(lw->flux_down, lw->device);
    }
    return GRTCODE_SUCCESS;
}


/*Calculate the spectral radiance [W*cm/m^2] using Planck's Law.*/
HOST DEVICE static int planck_law(fp_t const T, /*Temperature [K].*/
                                  fp_t const w, /*Wavenumber [cm-1].*/
                                  fp_t * const I /*Spectral radiance [W cm m-2].*/
                                 )
{
    /*Constants, defined as:
      c1 = 2*h*c*c with units: [W cm4 m-2]
      c2 = h*c/k with units: [cm K]
      where:
      h = planck's constant
      c = speed of light
      k = Boltzmann's constant*/
    fp_t const c1 = 1.1910429526245744e-8;
    fp_t const c2 = 1.4387773538277202;
    not_null(I);
    clear_floating_point_exceptions();
    fp_t e = c2*w/T;
    if (e > MAX_EXP_ARG)
    {
        /*Clamp down exponential argument to prevent overflow.*/
        e = MAX_EXP_ARG;
    }
    e = EXP(e);
    *I = (c1*w*w*w)/(e - 1.);
    catch_floating_point_exceptions();
    return GRTCODE_SUCCESS;
}


/*Calculate an effective spectral radiance [W cm m-2], taking into account a
  temperature difference between the center and edge of a layer.  This
  method was taken from https://doi.org/10.1029/92JD01419.*/
HOST DEVICE static int effective_planck(fp_t const Tcenter, /*Temperature [K] in the layer center.*/
                                        fp_t const Tedge, /*Temperature [K] at the layer edge.*/
                                        fp_t const w, /*Wavenumber [cm-1].*/
                                        fp_t const tau, /*Optical depth of layer at input wavenumber.*/
                                        fp_t * const I /*Spectral radiance [W cm m-2].*/
                                       )
{
    not_null(I);
    fp_t const a = 0.193; /*Using equation 16.*/
    fp_t const b = 0.013;
    fp_t planck_Tcenter;
    catch(planck_law(Tcenter, w, &planck_Tcenter));
    fp_t planck_Tedge;
    catch(planck_law(Tedge, w, &planck_Tedge));
    clear_floating_point_exceptions();
    *I = (planck_Tcenter + (a*tau + b*tau*tau)*planck_Tedge)/(1. + a*tau + b*tau*tau);
    catch_floating_point_exceptions();
    return GRTCODE_SUCCESS;
}


/*Calculate the upward and downward longwave radiative fluxes
  per wavenumber at each atmospheric pressure level.
  Arrays of size (0:n-1) must be stored so that index 0
  corresponds to top of the atmosphere, while index n-1 corresponds to
  the atmospheric layer/level nearest to the Earth's surface.  A four-
  stream approach is currently employed.*/
HOST DEVICE static int lw_flux(int const nlevels, /*Number of atmospheric pressure levels.*/
                               fp_t const w, /*Wavenumber [cm-1].*/
                               fp_t const T_surf, /*Temperature [K] of the Earth's surface.*/
                               fp_t const * const T_layers, /*Temperature [K] (layer).*/
                               fp_t const * const T_levels, /*Temperature [K] (level).*/
                               fp_t const * const tau, /*Optical depth (layer, wavenumber).*/
                               fp_t const emis, /*Surface emissivity (wavenumber).*/
                               fp_t * const flux_up, /*Upward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                               fp_t * const flux_down /*Downward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                              )
{
    num_levels_in_range(nlevels);
    wavenumber_in_range(w);
    temperature_in_range(T_surf);
    probability_in_range(emis);
    not_null(T_layers);
    not_null(T_levels);
    not_null(tau);
    not_null(flux_up);
    not_null(flux_down);
    clear_floating_point_exceptions();
    int const nlayers = nlevels - 1;
    int i;
    for (i=0; i<nlayers; ++i)
    {
        temperature_in_range(T_layers[i]);
        temperature_in_range(T_levels[i]);
        optical_depth_in_range(tau[i]);
    }
    temperature_in_range(T_levels[nlevels-1]);

    /*Set stream parameters.*/
    fp_t c1[4];
    c1[0] = -14.402613260847248;
    c1[1] = -3.0302159969901132;
    c1[2] = -1.4925584280108841;
    c1[3] = -1.0746123148178333;
    fp_t c2[4];
    c2[0] = 0.07587638482015649;
    c2[1] = 0.676114979733751;
    c2[2] = 1.3726594476601073;
    c2[3] = 1.0169418413757783;

    /*Zero out the flux arrays.*/
    memset(flux_down, 0, sizeof(*flux_down)*nlevels);
    memset(flux_up, 0, sizeof(*flux_up)*nlevels);

    /*Loop through the streams.*/
    int j;
    for (j=0; j<4; ++j)
    {
        /*Calculate the extinction for each layer.*/
        fp_t ext[MAX_NUM_LEVELS-1];
        for (i=0; i<nlayers; ++i)
        {
            fp_t e = c1[j]*tau[i];
            if (e > MAX_EXP_ARG)
            {
                /*Clamp down exponential argument to prevent overflow.*/
                e = MAX_EXP_ARG;
            }
            ext[i] = exp(e);
            catch_floating_point_exceptions();
        }

        /*Downward pass.*/
        fp_t I_down = 0.;
        for (i=0; i<nlayers; ++i)
        {
            fp_t val;
            catch(effective_planck(T_layers[i], T_levels[i+1], w, tau[i], &val));
            clear_floating_point_exceptions();
            fp_t p = (1. - ext[i])*val;
            I_down = p + I_down*ext[i];
            flux_down[i+1] += c2[j]*I_down;
            catch_floating_point_exceptions();
        }

        /*Upward pass.*/
        fp_t I_up;
        catch(planck_law(T_surf, w, &I_up));
        I_up = emis*I_up + (1 - emis)*I_down;
        flux_up[nlevels-1] += c2[j]*I_up;
        for (i=nlayers-1; i>=0; --i)
        {
            fp_t val;
            catch(effective_planck(T_layers[i], T_levels[i], w, tau[i], &val));
            clear_floating_point_exceptions();
            fp_t p = (1. - ext[i])*val;
            I_up = p + I_up*ext[i];
            flux_up[i] += c2[j]*I_up;
            catch_floating_point_exceptions();
        }
    }
    return GRTCODE_SUCCESS;
}


/*Calculate the longwave fluxes in a column at each wavenumber.*/
static int lw_fluxes_kernel(int const num_levels, /*Number of atmospheric pressure levels.*/
                            double const w0, /*Spectral grid lower bound [cm-1].*/
                            double const wres, /*Spectral grid resolution [cm-1].*/
                            uint64_t const num_wpoints, /*Spectral grid size.*/
                            fp_t const T_surf, /*Surface temperature [K].*/
                            fp_t const * const T_layers, /*Temperature [K] (layer).*/
                            fp_t const * const T_levels, /*Temperature [K] (level).*/
                            fp_t const * const tau, /*Optical depth (layer, wavenumber).*/
                            fp_t const * const omega, /*Single-scatter albedo (layer, wavenumber).*/
                            fp_t const * const emis, /*Surface emissivity (wavenumber).*/
                            fp_t * const flux_up, /*Upward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                            fp_t * const flux_down /*Downward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                           )
{
    uint64_t i;
#pragma omp parallel for schedule(static) default(shared) private(i)
    for (i=0; i<num_wpoints; ++i)
    {
        fp_t tau_buf[MAX_NUM_LEVELS];
        fp_t flux_up_buf[MAX_NUM_LEVELS];
        fp_t flux_down_buf[MAX_NUM_LEVELS];
        fp_t const w = w0 + i*wres;
        int const num_layers = num_levels - 1;
        int j;
        for (j=0; j<num_layers; ++j)
        {
            uint64_t o = j*num_wpoints+i;
            tau_buf[j] = tau[o]*(1. - omega[o]);
        }
        lw_flux(num_levels, w, T_surf, T_layers, T_levels, tau_buf, emis[i],
                flux_up_buf, flux_down_buf);
        for (j=0; j<num_levels; ++j)
        {
            flux_up[j*num_wpoints+i] = flux_up_buf[j];
            flux_down[j*num_wpoints+i] = flux_down_buf[j];
        }
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/*Calculate the longwave fluxes in a column at each wavenumber.*/
__global__ static void lw_fluxes_kernel_d(int const num_levels, /*Number of atmospheric pressure levels.*/
                                          double const w0, /*Spectral grid lower bound [cm-1].*/
                                          double const wres, /*Spectral grid resolution [cm-1].*/
                                          uint64_t const num_wpoints, /*Spectral grid size.*/
                                          fp_t const T_surf, /*Surface temperature [K].*/
                                          fp_t const * const T_layers, /*Temperature [K] (layer).*/
                                          fp_t const * const T_levels, /*Temperature [K] (level).*/
                                          fp_t const * const tau, /*Optical depth (layer, wavenumber).*/
                                          fp_t const * const omega, /* Single-scatter albedo (layer, wavenumber).*/
                                          fp_t const * const emis, /*Surface emissivity (wavenumber).*/
                                          fp_t * const flux_up, /*Upward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                                          fp_t * const flux_down /*Downward longwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                                         )
{
    uint64_t tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num_wpoints)
    {
        fp_t tau_buf[MAX_NUM_LEVELS];
        fp_t flux_up_buf[MAX_NUM_LEVELS];
        fp_t flux_down_buf[MAX_NUM_LEVELS];
        fp_t const w = w0 + tid*wres;
        int const num_layers = num_levels - 1;
        int j;
        for (j=0; j<num_layers; ++j)
        {
            uint64_t o = j*num_wpoints+tid;
            tau_buf[j] = tau[o]*(1. - omega[o]);
        }
        lw_flux(num_levels, w, T_surf, T_layers, T_levels, tau_buf, emis[tid],
                flux_up_buf, flux_down_buf);
        for (j=0; j<num_levels; ++j)
        {
            flux_up[j*num_wpoints+tid] = flux_up_buf[j];
            flux_down[j*num_wpoints+tid] = flux_down_buf[j];
        }
    }
    return;
}
#endif


/*Calculate the longwave radiative fluxes at each spectral grid point at each
  atmospheric level in the column.*/
EXTERN int calculate_lw_fluxes(Longwave_t * const lw, Optics_t const * const optics,
                               fp_t const T_surf, fp_t * const T_layers,
                               fp_t * const T_levels, fp_t * const emis,
                               fp_t * const flux_up, fp_t * const flux_down)
{
    not_null(lw);
    not_null(optics);
    not_null(T_layers);
    not_null(T_levels);
    not_null(emis);
    not_null(flux_up);
    not_null(flux_down);
    assert(lw->device, optics->device);
    assert(lw->num_levels, optics->num_layers+1);
    int same_grids;
    catch(compare_spectral_grids(&(lw->grid), &(optics->grid), &same_grids));
    assert(same_grids, 1);
    uint64_t n = lw->grid.n;
    if (lw->device == HOST_ONLY)
    {
        lw->layer_temperature = T_layers;
        lw->level_temperature = T_levels;
        lw->emissivity = emis;
        lw->flux_up = flux_up;
        lw->flux_down = flux_down;
    }
    else
    {
        gmemcpy(lw->layer_temperature, T_layers, optics->num_layers, lw->device, FROM_HOST);
        gmemcpy(lw->level_temperature, T_levels, lw->num_levels, lw->device, FROM_HOST);
        gmemcpy(lw->emissivity, emis, n, lw->device, FROM_HOST);
    }
    glaunch(lw_fluxes_kernel, n, lw->device, lw->num_levels, lw->grid.w0,
            lw->grid.dw, n, T_surf, lw->layer_temperature, lw->level_temperature,
            optics->tau, optics->omega, lw->emissivity, lw->flux_up, lw->flux_down);
    if (lw->device != HOST_ONLY)
    {
        gmemcpy(flux_up, lw->flux_up, n*lw->num_levels, lw->device, FROM_DEVICE);
        gmemcpy(flux_down, lw->flux_down, n*lw->num_levels, lw->device, FROM_DEVICE);
    }
    return GRTCODE_SUCCESS;
}
