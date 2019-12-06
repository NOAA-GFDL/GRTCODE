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

#ifndef CURTIS_GODSON_H_
#define CURTIS_GODSON_H_

#include "floating_point_type.h"


/** @brief Calculate integrated number densities.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_number_densities(int const num_layers, /**< Number of atmospheric layers.*/
                          fp_t const * const p, /**< Pressure [atm] (level).*/
                          fp_t * const n /**< Integrated number densities [cm-2] (layer).*/
                         );


#ifdef __NVCC__
/** @brief Calculate integrated number densities.*/
__global__ void calc_number_densities_d(int const num_layers, /**< Number of atmospheric layers.*/
                                        fp_t const * const p, /**< Pressure [atm] (level).*/
                                        fp_t * const n /**< Integrated number densities [cm-2] (layer).*/
                                       );
#endif


/** @brief Calculate layer pressures and temperatures.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_pressures_and_temperatures(int const num_layers, /**< Number of atmospheric layers.*/
                                    fp_t const * const p, /**< Pressure [atm] (level).*/
                                    fp_t const * const t, /**< Temperature [K] (level).*/
                                    fp_t * const pavg, /**< Pressure [atm] (layer).*/
                                    fp_t * const tavg /**< Pressure [atm] (layer).*/
                                   );


#ifdef __NVCC__
/** @brief Calculate layer pressures and temperatures.*/
__global__ void calc_pressures_and_temperatures_d(int const num_layers, /**< Number of atmospheric layers.*/
                                                  fp_t const * const p, /**< Pressure [atm] (level).*/
                                                  fp_t const * const t, /**< Temperature [K] (level).*/
                                                  fp_t * const pavg, /**< Pressure [atm] (layer).*/
                                                  fp_t * const tavg /**< Pressure [atm] (layer).*/
                                                 );
#endif


/** @brief Calculate partial pressures and number densities.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_partial_pressures_and_number_densities(int const num_layers, /**< Number of atmospheric layers.*/
                                                fp_t const * const p, /**< Pressure [atm] (level).*/
                                                fp_t const * const x, /**< Abundance (level).*/
                                                fp_t const * const n, /**< Integrated number densities [cm-2] (layer).*/
                                                fp_t * const ps, /**< Partial pressure [atm] (layer).*/
                                                fp_t * const ns /**< Integrated molecular number densities [cm-2] (layer).*/
                                               );


#ifdef __NVCC__
/** @brief Calculate partial pressures and number densities.*/
__global__ void calc_partial_pressures_and_number_densities_d(int const num_layers, /**< Number of atmospheric layers.*/
                                                              fp_t const * const p, /**< Pressure [atm] (level).*/
                                                              fp_t const * const x, /**< Abundance (level).*/
                                                              fp_t const * const n, /**< Integrated number densities [cm-2] (layer).*/
                                                              fp_t * const ps, /**< Partial pressure [atm] (layer).*/
                                                              fp_t * const ns /**< Integrated molecular number densities [cm-2] (layer).*/
                                                             );
#endif


#endif
