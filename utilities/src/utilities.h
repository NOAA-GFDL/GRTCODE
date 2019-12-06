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

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "extern.h"
#include "floating_point_type.h"


/** @brief Calculate the area under the curve between two points.*/
typedef fp_t(*Area1d_t)(fp_t const * const, fp_t const * const);


/** @brief Interpolation function between two points.*/
typedef int(*Sample1d_t)(fp_t const * const, fp_t const * const,
                         fp_t const * const, fp_t * const, size_t);

/** @brief Turn on bit in bit field.
    @return GRTCODE_SUCCESS or an error code.*/
int activate(uint64_t * const bit_field, /**< Bit field.*/
             int const index /**< Bit index (0 - 63).*/
            );


/** @brief Calculate an angstrom exponent.
    @return Angstrom exponent value.*/
fp_t angstrom_exponent(fp_t tau1, /**< Optical thickness at wavelength lambda1.*/
                       fp_t tau2, /**< Optical thickness at wavelength lambda2.*/
                       fp_t lambda1, /**< Wavelength.*/
                       fp_t lambda2 /**< Wavelength.*/
                      );


/** @brief Intepolate between two points using an "angstrom exponent".
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int angstrom_exponent_sample(fp_t const * const x, /**< Wavelengths to interpolate between (2).*/
                                    fp_t const * const y, /**< Optical thickness to interpolate between (2).*/
                                    fp_t const * const newx, /**< Wavelengths to interpolate to (n).*/
                                    fp_t * const newy, /**< Interpolated optical thicknesses (n).*/
                                    size_t n /**< Number of elements in newx and newy arrays.*/
                                   );


/** @brief Extrapolate using a constant value.
    @return GRTCODE_SUCCESS or an error code.*/
int constant_extrapolation(fp_t const * const x, /**< Unused.*/
                           fp_t const * const y, /**< Function value used for extrapolation.*/
                           fp_t const * const newx, /**< Unused.*/
                           fp_t * const newy, /**< Extrapolated function values (n).*/
                           size_t n /**< Number of elements in the newy array.*/
                          );


/** @brief Copy a string into a buffer, checking its length.
    @return GRTCODE_SUCCESS or an error code.*/
int copy_str(char * const dest, /**< Buffer the string is copied into.*/
             char const * const src, /**< String that is copied.*/
             size_t const len /**< Size of buffer.*/
            );


/** @brief Free malloced memory, with error checks.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int free_ptr(void ** const p /**< Address of pointer that is freed.*/
                   );


/** @brief Integrate piece-wise between adjacent values in an array.
    @return GRTCODE_SUCCESS or an error code.*/
int integrate2(fp_t const * const x, /**< Integration domain points (n).*/
               fp_t const * const y, /**< Function values that are integrated (n).*/
               size_t n, /**< Number of elements in x and y arrays.*/
               fp_t * const s, /**< Integral result.*/
               Area1d_t area /**< Function that calculates the area under the curve
                                  for adjacent points.*/
              );


/** @brief Interpolate piece-wise between adjacent values in an array.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int interpolate2(fp_t const * const x, /**< Domain points to interpolate from (n).*/
                        fp_t const * const y, /**< Function values to interpolate from (n).*/
                        size_t n, /**< Number of elements in the x and y arrays.*/
                        fp_t const * const newx, /**< Domain points to interpolate to (newn).*/
                        fp_t * const newy, /**< Interpolated function values (newn).*/
                        size_t newn, /**< Number of elements in the newx and newy arrays.*/
                        Sample1d_t interp, /**< Interpolation function between adjacent points.*/
                        Sample1d_t extrap /**< Extrapolation function for points outside the original domain.*/
                       );


/** @brief Check if bit is turned on in bit field.
    @return 0 if bit is off, else nonzero.*/
int is_active(uint64_t const bit_field, /**< Bit field.*/
              int const index /**< Bit index (0 - 63).*/
             );


/** @brief Linearly interpolate between two points.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int linear_sample(fp_t const * const x, /**< Domain points to interpolate from (2).*/
                         fp_t const * const y, /**< Function values to interpolate from (2).*/
                         fp_t const * const newx, /**< Domain points to interpolate to (n).*/
                         fp_t * const newy, /**< Interpolated function values (n).*/
                         size_t n /**< Number of elements in the newx and newy arrays.*/
                        );


/** @brief Malloc memory, with error checks.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int malloc_ptr(void ** const p, /**< Address of pointer that will hold the
                                            allocated memory.*/
                      size_t const num_bytes /**< Size of the allocation [bytes].*/
                     );


/** @brief Determine if the values in an array are monotonically increasing.
    @return 1 if monotonically increasing, else 0.*/
int monotonically_increasing(fp_t const * const x, /**< Data array (n).*/
                             size_t n /**< Number of elements in x array.*/
                            );


/** @brief Open a file, with error checks.
    @return GRTCODE_SUCCESS or an error code.*/
int open_file(FILE **file, /**< File handle.*/
              char const * const name, /**< File name.*/
              char const * const mode /**< File mode.*/
             );


/** @brief Convert a string to a double.
    @return GRTCODE_SUCCESS or an error code.*/
int to_double(char const * const s, /**< String value.*/
              double * const d /**< String converted to a double.*/
             );


/** @brief Convert a double to fp_t.
    @return GRTCODE_SUCCESS or an error code.*/
int to_fp_t(double const d, /**< Double value.*/
            fp_t * const f /**< Double value converted to a fp_t value.*/
           );


/** @brief Convert a string to an integer.
    @return GRTCODE_SUCCESS or an error code.*/
int to_int(char const * const s, /**< String value.*/
           int * const i /**< String converted to an integer.*/
          );


/** @brief Calculate the integral between two points using the trapezoid rule.
    @return Integral result.*/
fp_t trapezoid(fp_t const * const x, /**< Domain points (2).*/
               fp_t const * const y /**< Function values (2).*/
              );


#endif
