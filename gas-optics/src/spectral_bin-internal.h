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

#ifndef SPECTRAL_BIN_INTERNAL_H_
#define SPECTRAL_BIN_INTERNAL_H_

#include <stdint.h>
#include "grtcode_utilities.h"
#include "spectral_bin.h"


/** @brief Number of interpolation points per bin.  Here we mimic RFM and use a
     quadratic interpolation scheme, thus 3 points are required.*/
#define NIP 3


/** @brief Given spectral grid properties, define the spectral bins.
    @return GRTCODE_SUCCESS or an error code.*/
int create_spectral_bins(SpectralBins_t *bins, /**< Spectral bins.*/
                         int const num_layers, /**< Number of atmospheric layers.*/
                         double const w0, /**< Lower bound [cm-1] of spectral grid.*/
                         uint64_t const n, /**< Number of spectral grid points.*/
                         double const wres, /**< Resolution [cm-1] of spectral grid.*/
                         double const bin_width, /**< Desired width [cm-1] of each spectral bin.*/
                         Device_t const device /**< Device id.*/
                        );


/** @brief Free memory stored in the spectral bins.
    @return GRTCODE_SUCCESS or an error code.*/
int destroy_spectral_bins(SpectralBins_t *bins /**< Spectral bins.*/
                         );


#endif
