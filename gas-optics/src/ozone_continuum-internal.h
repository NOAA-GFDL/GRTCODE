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

#ifndef OZONE_CONTINUUM_INTERNAL_H_
#define OZONE_CONTINUUM_INTERNAL_H_

#include "grtcode_utilities.h"
#include "ozone_continuum.h"


/** @brief Read in the ozone continuum cross-sections.
    @return GRTCODE_SUCCESS or an error code.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc, /**< Continuum cross-sections.*/
                              char const * const o3_ctm_file, /**< Input data file.*/
                              SpectralGrid_t const grid, /**< Spectral grid.*/
                              Device_t const device /**< Device id.*/
                             );


/** @brief Free memory used for the ozone continuum cross-sections.
    @return GRTCODE_SUCCESS or an error code.*/
int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc /**< Continum cross-sections.*/
                              );


#endif
