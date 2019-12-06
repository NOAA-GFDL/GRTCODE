/** @file.*/
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

#ifndef SHORTWAVE_H_
#define SHORTWAVE_H_

#include <stdint.h>
#include "grtcode_utilities.h"


/** @brief Shortwave.*/
typedef struct Shortwave
{
    int num_levels; /**< Number of atmospheric levels.*/
    SpectralGrid_t grid; /**< Spectral grid.*/
    Device_t device; /**< Device.*/
    fp_t *solar_flux; /**< Incident solar flux [W cm m-2] (wavenumber).*/
    fp_t *sfc_alpha_dir; /**< Surface albedo for direct beam (wavenumber).*/
    fp_t *sfc_alpha_dif; /**< Surface albedo for diffuse beam (wavenumber).*/
    fp_t *flux_up; /**< Upward radiative flux [W cm m-2] (level, wavenumber).*/
    fp_t *flux_down; /**< Downward radiative flux [W cm m-2] (level, wavenumber).*/
} Shortwave_t;


/** @brief Reserve memory for the shortwave.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int create_shortwave(Shortwave_t * const sw, /**< Shortwave.*/
                            int const num_levels, /**< Number of atmospheric levels.*/
                            SpectralGrid_t const * const grid, /**< Spectral grid.*/
                            Device_t const * const device /**< Device.*/
                           );


/** @brief Free memory for the shortwave.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int destroy_shortwave(Shortwave_t * const sw /**< Shortwave.*/
                            );


/** @brief Calculate the shortwave radiative fluxes at each spectral grid point at
           each atmospheric level in the column.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int calculate_sw_fluxes(Shortwave_t * const sw, /**< Shortwave object.*/
                               Optics_t const * const optics, /**< Optics object.*/
                               fp_t const mu_dir, /**< Cosine of zenith angle for direct beam.*/
                               fp_t const mu_dif, /**< Cosine of zenith angle for diffuse beam.*/
                               fp_t * const sfc_alpha_dir, /**< Surface albedo for direct beam (wavenumber).*/
                               fp_t * const sfc_alpha_dif, /**< Surface albedo for diffuse beam (wavenumber).*/
                               fp_t const total_solar_irradiance, /**< Total solar irradiance [W m-2].*/
                               fp_t * const solar_flux, /**< Solar flux [cm] (wavenumber).*/
                               fp_t * const flux_up, /**< Upward flux [W cm m-2] (level, wavenumber).*/
                               fp_t * const flux_down /**< Downward flux [W cm m-2] (level, wavenumber).*/
                              );


#endif
