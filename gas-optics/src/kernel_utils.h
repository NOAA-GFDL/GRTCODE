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

#ifndef KERNEL_UTILS_H_
#define KERNEL_UTILS_H_

#include <stdint.h>
#include "debug.h"
#include "grtcode_utilities.h"


/** @brief Get the indices of an array that bracket a value.  The input
           array must be sorted.
    @return GRTCODE_SUCCESS or an error code.*/
HOST DEVICE int bracket(uint64_t const array_size, /**< Size of input array.*/
                        fp_t const * const array, /**< Input array.*/
                        fp_t const val, /**< Value to bracket.*/
                        uint64_t * const left, /**< Index of maximum array value < val.*/
                        uint64_t * const right /**< Index of minimum array value > val.*/
                       );


/** @brief Perform a quadratic interpolation, but set all values less than
           zero equal to zero.
    @return GRTCODE_SUCCESS or an error code.*/
HOST DEVICE int bin_quad_interp(fp_t const * const x, /**< Domain points to interpolate from.*/
                                fp_t const * const y, /**< Function values to interpolate from.*/
                                uint64_t const left, /**< Starting index.*/
                                uint64_t const right, /**< Ending index.*/
                                fp_t const w0, /**< Spectral grid lower bound [cm-1].*/
                                fp_t const dw, /**< Spectral grid spacing [cm-1].*/
                                fp_t * const tau /**< Interpolated function values.*/
                               );


/** @brief Copy optical depths from the coarse mesh to the fine mesh in a bin.
    @return GRTCODE_SUCCESS or an error code.*/
HOST DEVICE int bin_no_interp(uint64_t const left, /**< Starting index.*/
                              uint64_t const right, /**< Ending index.*/
                              fp_t const * const taub, /**< Data to be copied.*/
                              fp_t * const tau /**< Buffer to copy into.*/
                             );


#endif
