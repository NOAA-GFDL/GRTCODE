/** @file*/
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

#ifndef SOLAR_FLUX_H_
#define SOLAR_FLUX_H_

#include "grtcode_utilities.h"


/** @brief Solar flux.*/
typedef struct SolarFlux
{
    SpectralGrid_t grid; /**< Spectral grid.*/
    fp_t *incident_flux; /**< Incident solar flux [cm] (wavenumber).*/
    uint64_t n; /**< Size of spectral grid.*/
} SolarFlux_t;


/** @brief Read in data for the solar flux.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int create_solar_flux(SolarFlux_t * const solar_flux, /**< Solar flux.*/
                             SpectralGrid_t const * const grid, /**< Spectral grid.*/
                             char const * const filepath /**< Solar flux csv file.*/
                            );


/** @brief Free memory for the solar flux.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int destroy_solar_flux(SolarFlux_t * const solar_flux /**< Solar flux.*/
                             );


#endif
