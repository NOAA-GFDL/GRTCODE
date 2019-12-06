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

#ifndef GRTCODE_CONFIG_H_
#define GRTCODE_CONFIG_H_

#include <float.h>


#ifdef SINGLE_PRECISION
#define ABS fabsf
#define EXP expf
#define LOG logf
/*ln(max(float)) = 88.72284.  Let's use 80 so we have some runway.*/
#define MAX_EXP_ARG 80.f
#define MAX_OPTICAL_DEPTH FLT_MAX
#define MAX_SOLAR_FLUX FLT_MAX
#define MIN_COS_ZENITH_ANGLE FLT_MIN
#define POW powf
#define SQRT sqrtf
#else
#define ABS fabs
#define EXP exp
#define LOG log
/*ln(max(double)) = 709.782712893384.  Let's use 700 so we have some runway.*/
#define MAX_EXP_ARG 700.
#define MAX_OPTICAL_DEPTH DBL_MAX
#define MAX_SOLAR_FLUX DBL_MAX
#define MIN_COS_ZENITH_ANGLE DBL_MIN
#define POW pow
#define SQRT sqrt
#endif


#define DEFAULT_GPU 0
#define MIN_NUM_LAYERS 1
#define MAX_NUM_LAYERS 200
#define MIN_NUM_LEVELS MIN_NUM_LAYERS+1
#define MAX_NUM_LEVELS MAX_NUM_LAYERS+1
#define MIN_WAVENUMBER 1.
#define MAX_WAVENUMBER 50000.
#define MIN_RESOLUTION 0.001
#define MAX_RESOLUTION 10.
#define MIN_PROBABILITY 0.
#define MAX_PROBABILITY 1.
#define MIN_TEMPERATURE 100.
#define MAX_TEMPERATURE 500.
#define MIN_OPTICAL_DEPTH 0.
#define MAX_COS_ZENITH_ANGLE 1.
#define MIN_SOLAR_FLUX 0.
#define MIN_ASYMMETRIC_FACTOR -1.
#define MAX_ASYMMETRIC_FACTOR 1.


#define probability_in_range(p) \
    in_range(p, MIN_PROBABILITY, MAX_PROBABILITY);


#define num_levels_in_range(n) \
    in_range(n, MIN_NUM_LEVELS, MAX_NUM_LEVELS);


#define wavenumber_in_range(w) \
    in_range(w, MIN_WAVENUMBER, MAX_WAVENUMBER);


#define temperature_in_range(T) \
    in_range(T, MIN_TEMPERATURE, MAX_TEMPERATURE);


#define optical_depth_in_range(t) \
    in_range(t, MIN_OPTICAL_DEPTH, MAX_OPTICAL_DEPTH);


#define cosine_zenith_angle_in_range(c) \
    in_range(c, MIN_COS_ZENITH_ANGLE, MAX_COS_ZENITH_ANGLE);


#define solar_flux_in_range(f) \
    in_range(f, MIN_SOLAR_FLUX, MAX_SOLAR_FLUX);


#define asymmetric_factor_in_range(g) \
    in_range(g, MIN_ASYMMETRIC_FACTOR, MAX_ASYMMETRIC_FACTOR);


#endif
