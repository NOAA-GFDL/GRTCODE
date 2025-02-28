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
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "shortwave.h"


/*Reserve memory for the shortwave.*/
EXTERN int create_shortwave(Shortwave_t * const sw, int const num_levels,
                            SpectralGrid_t const * const grid, Device_t const * const device)
{
    not_null(sw);
    not_null(grid);
    not_null(device);
    in_range(num_levels, MIN_NUM_LEVELS, MAX_NUM_LEVELS);
    sw->num_levels = num_levels;
    sw->grid = *grid;
    sw->device = *device;
    if (sw->device != HOST_ONLY)
    {
        gmalloc(sw->solar_flux, sw->grid.n, sw->device);
        gmalloc(sw->sfc_alpha_dir, sw->grid.n, sw->device);
        gmalloc(sw->sfc_alpha_dif, sw->grid.n, sw->device);
        gmalloc(sw->flux_up, sw->grid.n*num_levels, sw->device);
        gmalloc(sw->flux_down, sw->grid.n*num_levels, sw->device);
    }
    return GRTCODE_SUCCESS;
}


/*Free memory for the shortwave.*/
EXTERN int destroy_shortwave(Shortwave_t * const sw)
{
    not_null(sw);
    if (sw->device != HOST_ONLY)
    {
        gfree(sw->solar_flux, sw->device);
        gfree(sw->sfc_alpha_dir, sw->device);
        gfree(sw->sfc_alpha_dif, sw->device);
        gfree(sw->flux_up, sw->device);
        gfree(sw->flux_down, sw->device);
    }
    return GRTCODE_SUCCESS;
}


/*Perform scaling required by the delta_Eddington Method.  The scaling parameters are
  described in https://doi.org/10.1175/1520-0469(1976)033<2452:TDEAFR>2.0.CO;2.*/
HOST DEVICE static int delta_eddington_scaling_jww1976(fp_t const omega, /*Single-scatter albedo.*/
                                                       fp_t const g, /*Asymmetry factor.*/
                                                       fp_t const tau, /*Optical depth.*/
                                                       fp_t * const omega_s, /*Scaled single-scatter albedo.*/
                                                       fp_t * const g_s, /*Scaled asymmetry factor.*/
                                                       fp_t * const tau_s /*Scaled optical depth.*/
                                                      )
{
    /*Check inputs.*/
    probability_in_range(omega);
    asymmetric_factor_in_range(g);
    optical_depth_in_range(tau);
    not_null(omega_s);
    not_null(g_s);
    not_null(tau_s);

    /*Calculate scaled quantities, divide-by-zero for g = -1.*/
    clear_floating_point_exceptions();
    *g_s = g/(g + 1.); /*Equation 5b.*/
    fp_t const f = g*g; /*Equation 5a.*/
    *omega_s = (1.-f)*omega/(1.-omega*f); /*Equation 14.*/
    *tau_s = tau*(1.-omega*f); /*Equation 13.*/
    catch_floating_point_exceptions();
    return GRTCODE_SUCCESS;
}


/*Calcluate the reflectivity and transmittance of a layer using the method described in
  https://doi.org/10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2.*/
HOST DEVICE static int meador_weaver_1980(fp_t const omega, /*Scaled single-scatter albedo.*/
                                          fp_t const tau, /*Optical depth.*/
                                          fp_t const mu, /*Cosine of the beam's zenith angle.*/
                                          fp_t const gamma1, /*See paper.*/
                                          fp_t const gamma2, /*See paper.*/
                                          fp_t const gamma3, /*See paper.*/
                                          fp_t * const R, /*Layer reflectivity.*/
                                          fp_t * const T, /*Layer transmittance.*/
                                          fp_t * const T_pure /*Amount of incident radiation
                                                                that passes through the layer
                                                                without being absorbed or scattered.*/
                                         )
{
    not_null(R);
    not_null(T);
    clear_floating_point_exceptions();

    if (omega <= 0.0)
    {
        /*There is no scattering, thus no reflection.*/
        *R = 0.;
        *T = exp(-tau/mu);
        if (T_pure != NULL)
        {
            *T_pure = *T;
        }
    }
    else
    {
        /*Calculate more parameters defined in the paper.*/
        double const gamma4 = 1. - gamma3; /*Equation 21.*/
        double const alpha1 = gamma1*gamma4 + gamma2*gamma3; /*Equation 16.*/
        double const alpha2 = gamma1*gamma3 + gamma2*gamma4; /*Equation 17.*/
        double const k = sqrt(gamma1*gamma1 - gamma2*gamma2); /*Equation 18.*/
        catch_floating_point_exceptions();

        /*Clamp down optical depth values to prevent numerical issues when
          the optical depth is very large.  This can lead to some relatively
          large errors in the transmittance (although they are large errors
          in very small values).*/
        fp_t t = tau;
        if (1./mu > k && tau/mu > MAX_EXP_ARG)
        {
            t = MAX_EXP_ARG*mu;
        }
        else if (tau*k > MAX_EXP_ARG)
        {
            t = MAX_EXP_ARG/k;
        }

        fp_t const tp = exp(t/mu);
        catch_floating_point_exceptions();
        if (tp <= 1.0)
        {
            /*There is no gas in the layer.*/
            *R = 0.;
            *T = 1.;
            if (T_pure != NULL)
            {
                *T_pure = 1.;
            }
        }
        else
        {
            fp_t const tm = exp(-t/mu);
            fp_t const tkm = exp(-t*k);
            fp_t const tkp = exp(t*k);
            catch_floating_point_exceptions();
            if (T_pure != NULL)
            {
                *T_pure = tm;
            }
            if (omega >= 1.)
            {
                /*Scattering is conservative (i.e., there is no absorption.).*/
                *R = (1./(1.+gamma1*t))*(gamma1*t + (gamma3 - gamma1*mu)*
                     (1.-tm)); /*Equation 24.*/
                catch_floating_point_exceptions();
                probability_in_range(*R);
                *T = 1. - *R; /*Equation 24.*/
                catch_floating_point_exceptions();
                probability_in_range(*T);
            }
            else
            {
                /*Equation 14.*/
                *R = (omega/((1.-k*k*mu*mu)*((k+gamma1)*tkp + (k-gamma1)*tkm)))*
                     ((1.-k*mu)*(alpha2+k*gamma3)*tkp - (1.+k*mu)*(alpha2-k*gamma3)*tkm -
                     2.*k*(gamma3-alpha2*mu)*tm);
                catch_floating_point_exceptions();
                probability_in_range(*R);

                /*Equation 15.*/
                *T = tm*(1. - (omega/((1.-k*k*mu*mu)*((k+gamma1)*tkp +
                     (k-gamma1)*tkm)))*((1.+k*mu)*(alpha1+k*gamma4)*tkp -
                     (1.-k*mu)*(alpha1-k*gamma4)*tkm - 2.*k*(gamma4 + alpha1*mu)*tp));
                catch_floating_point_exceptions();
                probability_in_range(*T);
            }
        }
    }

    if (T_pure != NULL)
    {
        if (*T_pure > *T)
        {
            *T = *T_pure;
        }
    }
    return GRTCODE_SUCCESS;
}


/*Use the Eddington approximation to calculate the reflectivity and
  transmittance of a plane-parallel atmospheric layer.  The specific
  form of the Eddinton approximation used here is described in
  https://doi.org/10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2.*/
HOST DEVICE static int eddington_mw1980(fp_t const omega, /*Scaled single-scatter albedo.*/
                                        fp_t const tau, /*Scaled layer optical depth.*/
                                        fp_t const mu, /*Cosine of the solar zenith angle.*/
                                        fp_t const g, /*Scaled asymmetry factor.*/
                                        fp_t * const R, /*Layer reflectivity.*/
                                        fp_t * const T, /*Layer transmittance.*/
                                        fp_t * const T_pure /*Amount of incident radiation that
                                                              passes through the layer without
                                                              being absorbed or scattered.*/
                                       )
{
    /*Calculate the Eddington approximation parameters, as defined in
      row 1 of table 1.*/
    clear_floating_point_exceptions();
    fp_t const gamma1 = 0.25*(7. - omega*(4. + 3.*g));
    fp_t const gamma2 = -0.25*(1. - omega*(4. - 3.*g));
    fp_t const gamma3 = 0.25*(2. - 3.*g*mu);
    catch_floating_point_exceptions();

    /*Calculate the reflectivity and transmittance of the layer.*/
    catch(meador_weaver_1980(omega, tau, mu, gamma1, gamma2, gamma3, R, T, T_pure));
    return GRTCODE_SUCCESS;
}


/*Calculate the shortwave upward and downward fluxes using the "adding"
  method.  An overview of this method is provided in the appendix of
  https://doi.org/10.1029/94JD01310.*/
HOST DEVICE static int sw_adding(int const nlevels, /*Number of atmospheric pressure levels.*/
                                 fp_t const * const R_dir, /*Reflectivity of the direct beam (layer).*/
                                 fp_t const * const R_dif, /*Reflectivity of the diffuse beam (layer).*/
                                 fp_t const * const T_dir, /*Transmittance of the direct beam (layer).*/
                                 fp_t const * const T_dif, /*Transmittance of the diffuse beam (layer).*/
                                 fp_t const * const T_pure, /*Amount of incident radiation that passes
                                                              through each layer without being absorbed
                                                              or scattered.*/
                                 fp_t const sfc_alpha_dir, /*Albedo of the Earth's surface for a direct beam.*/
                                 fp_t const sfc_alpha_dif, /*Albedo of the Earth's surface for a diffuse beam.*/
                                 fp_t * const R, /*Total upward reflectance (level).*/
                                 fp_t * const T /*Total downward reflectance (level).*/
                                )
{
    not_null(R_dir);
    not_null(R_dif);
    not_null(T_dir);
    not_null(T_dif);
    not_null(T_pure);
    not_null(R);
    not_null(T);

    int const nlayers = nlevels - 1;
    fp_t R_dir_downward[MAX_NUM_LEVELS]; /*Reflectance for a downward traveling direct beam (layer).*/
    fp_t R_dif_downward[MAX_NUM_LEVELS]; /*Reflectance for a downward traveling diffuse beam (layer).*/
    fp_t R_dif_upward[MAX_NUM_LEVELS]; /*Reflectance for an upward traveling diffuse beam (layer).*/

    /*For a downward traveling beam incident on a layer from above,
      calculate the reflectance of the top of the layer.  Start at the
      lowest layer in the atmosphere, and then build upward.*/
    R_dir_downward[nlevels-1] = sfc_alpha_dir;
    R_dif_downward[nlevels-1] = sfc_alpha_dif;
    int i;
    for (i=nlayers-1; i>=0; --i)
    {
        clear_floating_point_exceptions();
        fp_t A = T_pure[i];
        fp_t B = 1./(1. - R_dif[i]*R_dif_downward[i+1]);
        R_dir_downward[i] = R_dir[i] + (A*R_dir_downward[i+1] +
                            (T_dir[i] - A)*R_dif_downward[i+1])*T_dif[i]*B;
        catch_floating_point_exceptions();
        probability_in_range(R_dir_downward[i]);

        clear_floating_point_exceptions();
        R_dif_downward[i] = R_dif[i] + T_dif[i]*T_dif[i]*R_dif_downward[i+1]*B;
        catch_floating_point_exceptions();
        probability_in_range(R_dif_downward[i]);
    }

    /*For an upward traveling beam incident on a layer from below,
      calculate the reflectance of the bottom of the layer.  Start at the
      top layer in the atmosphere, and then build downward.*/
    R_dif_upward[0] = R_dif[0];
    for (i=1; i<nlayers; ++i)
    {
        clear_floating_point_exceptions();
        fp_t B = 1./(1. - R_dif[i]*R_dif_upward[i-1]);
        R_dif_upward[i] = R_dif[i] + T_dif[i]*T_dif[i]*R_dif_upward[i-1]*B;
        catch_floating_point_exceptions();
        probability_in_range(R_dif_upward[i]);
    }

    /*Calculate the flux through each pressure level.*/
    fp_t solar_flux = 1.;
    fp_t dir_beam = solar_flux;
    fp_t dif_beam = 0.;
    R[0] = dir_beam*R_dir_downward[0];
    T[0] = dir_beam;
    for (i=1; i<nlevels; ++i)
    {
        if (i > 1)
        {
            fp_t C = 1./(1.-R_dif[i-1]*R_dif_upward[i-2]);
            dif_beam = (dir_beam*R_dir[i-1]*R_dif_upward[i-2] + dif_beam)*T_dif[i-1]*C +
                       dir_beam*(T_dir[i-1]-T_pure[i-1]);
        }
        else
        {
            dif_beam = dir_beam*(T_dir[i-1]-T_pure[i-1]);
        }
        in_range(dif_beam, 0., solar_flux);
        dir_beam *= T_pure[i-1];
        in_range(dir_beam, 0., solar_flux);
        fp_t B = 1./(1. - R_dif_downward[i]*R_dif_upward[i-1]);
        R[i] = (dir_beam*R_dir_downward[i] + dif_beam*R_dif_downward[i])*B;
        T[i] = dir_beam*(1. + R_dir_downward[i]*R_dif_upward[i-1]*B) + dif_beam*B;
    }
    return GRTCODE_SUCCESS;
}


/*Calculate the upward and downward shortwave radiative fluxes
  per wavenumber at each atmospheric pressure level.
  Arrays of size (0:n-1) must be stored so that index 0
  corresponds to top of the atmosphere, while index n-1 corresponds to
  the atmospheric layer/level nearest to the Earth's surface.  This
  routine currently employs the delta-Eddington and adding methods.*/
HOST DEVICE static int sw_flux(int const nlevels, /*Number of atmospheric pressure levels.*/
                               fp_t const * const omega, /*Single-scatter albedo (layer).*/
                               fp_t const * const g, /*Asymmetry factor (layer).*/
                               fp_t const * const tau, /*Optical depth (layer).*/
                               fp_t const mu_dir, /*Cosine of the zenith angle for a direct beam.*/
                               fp_t const mu_dif, /*Cosine of the zenith angle for a diffuse beam.*/
                               fp_t const sfc_alpha_dir, /*Albedo of the Earth's surface for a direct beam.*/
                               fp_t const sfc_alpha_dif, /*Albedo of the Earth's surface for a diffuse beam.*/
                               fp_t const solar_flux, /*Incident solar flux [cm] at the top of the atmosphere.*/
                               fp_t * const flux_up, /*Upward shortwave radiative flux [cm] (level).*/
                               fp_t * const flux_down /*Downward shortwave radiative flux [cm] (level).*/
                              )
{
    /*Check inputs.*/
    num_levels_in_range(nlevels);
    cosine_zenith_angle_in_range(mu_dir);
    cosine_zenith_angle_in_range(mu_dif);
    probability_in_range(sfc_alpha_dir);
    probability_in_range(sfc_alpha_dif);
    not_null(omega);
    not_null(g);
    not_null(tau);
    not_null(flux_up);
    not_null(flux_down);

    /*Local buffers.*/
    int const nlayers = nlevels - 1;
    fp_t R_dir[MAX_NUM_LEVELS-1]; /*Reflectivity for a direct beam (layer).*/
    fp_t R_dif[MAX_NUM_LEVELS-1]; /*Reflectivity for a diffuse beam (layer).*/
    fp_t T_dir[MAX_NUM_LEVELS-1]; /*Transmittance for a direct beam (layer).*/
    fp_t T_dif[MAX_NUM_LEVELS-1]; /*Transmittance for a diffuse beam (layer).*/
    fp_t T[MAX_NUM_LEVELS-1]; /*Amount of the incident beam that makes it through
                                without being absorbed or scattered (layer).*/

    int i;
    for (i=0; i<nlayers; ++i)
    {
        /*Perform the delta-Eddington scaling.*/
        fp_t omega_scaled;
        fp_t g_scaled;
        fp_t tau_scaled;
        catch(delta_eddington_scaling_jww1976(omega[i], g[i], tau[i], &omega_scaled,
                                              &g_scaled, &tau_scaled));

        /*Calculate the layer reflectivities and transmittances using the
          Eddington approximation for the direct beam.*/
        catch(eddington_mw1980(omega_scaled, tau_scaled, mu_dir, g_scaled,
                               &(R_dir[i]), &(T_dir[i]), &(T[i])));

        /*Calculate the layer reflectivities and transmittances using the
          Eddington approximation for the diffuse beam.*/
        catch(eddington_mw1980(omega_scaled, tau_scaled, mu_dif, g_scaled,
                               &(R_dif[i]), &(T_dif[i]), NULL));
    }

    /*Calculate the total reflectance and transmittance at each pressure
      level.*/
    catch(sw_adding(nlevels, R_dir, R_dif, T_dir, T_dif, T, sfc_alpha_dir,
                    sfc_alpha_dif, flux_up, flux_down));

    /*Calculate the upward and downward fluxes.*/
    for (i=0; i<nlevels; ++i)
    {
        flux_up[i] *= solar_flux*mu_dir;
        flux_down[i] *= solar_flux*mu_dir;
    }
    return GRTCODE_SUCCESS;
}


/*Calculate the shortwave fluxes in a column at each wavenumber.*/
static int sw_fluxes_kernel(int const num_levels, /*Number of atmospheric pressure levels.*/
                            uint64_t const num_wpoints, /*Spectral grid size.*/
                            fp_t const * const omega, /*Single-scatter albedo (layer, wavenumber).*/
                            fp_t const * const g, /*Asymmetry factor (layer, wavenumber).*/
                            fp_t const * const tau, /*Optical depth (layer, wavenumber).*/
                            fp_t const mu_dir, /*Cosine of the zenith angle for a direct beam.*/
                            fp_t const mu_dif, /*Cosine of the zenith angle for a diffuse beam.*/
                            fp_t const * const sfc_alpha_dir, /*Albedo of the Earth's surface for a direct beam (wavenumber).*/
                            fp_t const * const sfc_alpha_dif, /*Albedo of the Earth's surface for a diffuse beam (wavenumber).*/
                            fp_t const total_solar_irradiance, /*Total solar irradiance [W m-2].*/
                            fp_t const * const solar_flux, /*Incident solar flux [cm] at the top of the atmosphere (wavenumber).*/
                            fp_t * const flux_up, /*Upward shortwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                            fp_t * const flux_down /*Downward shortwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                           )
{
    uint64_t i;
#pragma omp parallel for schedule(static) default(shared) private(i)
    for (i=0; i<num_wpoints; ++i)
    {
        fp_t omega_buf[MAX_NUM_LEVELS];
        fp_t g_buf[MAX_NUM_LEVELS];
        fp_t tau_buf[MAX_NUM_LEVELS];
        fp_t flux_up_buf[MAX_NUM_LEVELS];
        fp_t flux_down_buf[MAX_NUM_LEVELS];
        int num_layers = num_levels - 1;
        int j;
        for (j=0; j<num_layers; ++j)
        {
            uint64_t offset = j*num_wpoints + i;
            omega_buf[j] = omega[offset];
            g_buf[j] = g[offset];
            tau_buf[j] = tau[offset];
        }
        sw_flux(num_levels, omega_buf, g_buf, tau_buf, mu_dir, mu_dif, sfc_alpha_dir[i],
                sfc_alpha_dif[i], solar_flux[i], flux_up_buf, flux_down_buf);
        for (j=0; j<num_levels; ++j)
        {
            uint64_t offset = j*num_wpoints + i;
            flux_up[offset] = total_solar_irradiance*flux_up_buf[j];
            flux_down[offset] = total_solar_irradiance*flux_down_buf[j];
        }
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/*Calculate the shortwave fluxes in a column at each wavenumber.*/
__global__ static void sw_fluxes_kernel_d(int const num_levels, /*Number of atmospheric pressure levels.*/
                                          uint64_t const num_wpoints, /*Spectral grid size.*/
                                          fp_t const * const omega, /*Single-scatter albedo (layer, wavenumber).*/
                                          fp_t const * const g, /*Asymmetry factor (layer, wavenumber).*/
                                          fp_t const * const tau, /*Optical depth (layer, wavenumber).*/
                                          fp_t const mu_dir, /*Cosine of the zenith angle for a direct beam.*/
                                          fp_t const mu_dif, /*Cosine of the zenith angle for a diffuse beam.*/
                                          fp_t const * const sfc_alpha_dir, /*Albedo of the Earth's surface for a direct beam (wavenumber).*/
                                          fp_t const * const sfc_alpha_dif, /*Albedo of the Earth's surface for a diffuse beam (wavenumber).*/
                                          fp_t const total_solar_irradiance, /*Total solar irradiance [W m-2].*/
                                          fp_t const * const solar_flux, /*Incident solar flux [cm] at the top of the atmosphere (wavenumber).*/
                                          fp_t * const flux_up, /*Upward shortwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                                          fp_t * const flux_down /*Downward shortwave radiative fluxes [W cm m-2] (level, wavenumber).*/
                                         )
{
    uint64_t tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < num_wpoints)
    {
        fp_t omega_buf[MAX_NUM_LEVELS];
        fp_t g_buf[MAX_NUM_LEVELS];
        fp_t tau_buf[MAX_NUM_LEVELS];
        fp_t flux_up_buf[MAX_NUM_LEVELS];
        fp_t flux_down_buf[MAX_NUM_LEVELS];
        int num_layers = num_levels - 1;
        int j;
        for (j=0; j<num_layers; ++j)
        {
            uint64_t offset = j*num_wpoints + tid;
            omega_buf[j] = omega[offset];
            g_buf[j] = g[offset];
            tau_buf[j] = tau[offset];
        }
        sw_flux(num_levels, omega_buf, g_buf, tau_buf, mu_dir, mu_dif, sfc_alpha_dir[tid],
                sfc_alpha_dif[tid], solar_flux[tid], flux_up_buf, flux_down_buf);
        for (j=0; j<num_levels; ++j)
        {
            uint64_t offset = j*num_wpoints + tid;
            flux_up[offset] = total_solar_irradiance*flux_up_buf[j];
            flux_down[offset] = total_solar_irradiance*flux_down_buf[j];
        }
    }
    return;
}
#endif


/*Calculate the shortwave radiative fluxes at each spectral grid point at
  each atmospheric level in the column.*/
EXTERN int calculate_sw_fluxes(Shortwave_t * const sw, Optics_t const * const optics,
                               fp_t const mu_dir, fp_t const mu_dif,
                               fp_t * const sfc_alpha_dir, fp_t * const sfc_alpha_dif,
                               fp_t const total_solar_irradiance,
                               fp_t * const solar_flux, fp_t * const flux_up,
                               fp_t * const flux_down)
{
    not_null(sw);
    not_null(optics);
    not_null(solar_flux);
    not_null(flux_up);
    not_null(flux_down);
    assert(sw->device, optics->device);
    assert(sw->num_levels, optics->num_layers+1);
    int same_grids;
    catch(compare_spectral_grids(&(sw->grid), &(optics->grid), &same_grids));
    assert(same_grids, 1);
    if (sw->device == HOST_ONLY)
    {
        sw->solar_flux = solar_flux;
        sw->sfc_alpha_dir = sfc_alpha_dir;
        sw->sfc_alpha_dif = sfc_alpha_dif;
        sw->flux_up = flux_up;
        sw->flux_down = flux_down;
    }
    else
    {
        gmemcpy(sw->solar_flux, solar_flux, sw->grid.n, sw->device, FROM_HOST);
        gmemcpy(sw->sfc_alpha_dir, sfc_alpha_dir, sw->grid.n, sw->device, FROM_HOST);
        gmemcpy(sw->sfc_alpha_dif, sfc_alpha_dif, sw->grid.n, sw->device, FROM_HOST);
    }
    glaunch(sw_fluxes_kernel, sw->grid.n, sw->device, sw->num_levels, sw->grid.n,
            optics->omega, optics->g, optics->tau, mu_dir, mu_dif, sw->sfc_alpha_dir,
            sw->sfc_alpha_dif, total_solar_irradiance, sw->solar_flux, sw->flux_up,
            sw->flux_down);
    if (sw->device != HOST_ONLY)
    {
        gmemcpy(flux_up, sw->flux_up, sw->grid.n*sw->num_levels, sw->device, FROM_DEVICE);
        gmemcpy(flux_down, sw->flux_down, sw->grid.n*sw->num_levels, sw->device, FROM_DEVICE);
    }
    return GRTCODE_SUCCESS;
}
