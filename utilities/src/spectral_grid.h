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

#ifndef SPECTRAL_GRID_H_
#define SPECTRAL_GRID_H_

#include <stddef.h>
#include <stdint.h>
#include "device.h"
#include "extern.h"
#include "floating_point_type.h"
#include "utilities.h"


/** @brief Spectral grid properties.*/
typedef struct SpectralGrid
{
    double dw; /**< Grid spacing [cm-1].*/
    uint64_t n; /**< Number of grid points.*/
    double wn; /**< Upper bound [cm-1].*/
    double w0; /**< Lower bound [cm-1].*/
} SpectralGrid_t;


/** @brief Determine if two spectral grids are the same.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int compare_spectral_grids(SpectralGrid_t const * const one, /**< Spectral grid.*/
                                  SpectralGrid_t const * const two, /**< Spectral grid.*/
                                  int * const result /**< 1 if they match, 0 if not.*/
                                 );


/** @brief Initialize a spectral grid.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int create_spectral_grid(SpectralGrid_t * const grid, /**< Spectral grid.*/
                                double const w0, /**< Lower bound [cm-1].*/
                                double const wn, /**< Upper bound [cm-1].*/
                                double const dw /**< Grid spacing [cm-1].*/
                               );


/** @brief Get the index of a grid point.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int grid_point_index(SpectralGrid_t const grid, /**< Spectral grid.*/
                            double const w, /**< Point to search for [cm-1].*/
                            uint64_t * const index /**< Grid index.*/
                           );


/** @brief Create an array of all grid points.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int grid_points(SpectralGrid_t const grid, /**< Spectral grid.*/
                       fp_t **buffer, /**< Grid points buffer [cm-1] (wavenumber).*/
                       Device_t const device /**< Device id.*/
                      );


/** @brief Interpolate an array onto the input spectral grid.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int interpolate_to_grid(SpectralGrid_t const grid, /**< Spectral grid.*/
                               fp_t const * const x, /**< Domain points to interpolate from [cm-1] (n).*/
                               fp_t const * const y, /**< Function values to interpolate from (n).*/
                               size_t const n, /**< Number of elements in the x and y arrays.*/
                               fp_t * const newy, /**< Interpolated function values (wavenumber).*/
                               Sample1d_t interp, /**< Interpolation function between adjacent points.*/
                               Sample1d_t extrap /**< Extrapolation function for points outside the original domain.*/
                              );


#endif
