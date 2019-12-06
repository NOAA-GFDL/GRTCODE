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

#ifndef CLOUD_OPTICS_H_
#define CLOUD_OPTICS_H_

#include "grtcode_utilities.h"


/** @brief Liquid cloud parameterization.*/
typedef int(*LiquidCloud_t)(fp_t const, fp_t const, SpectralGrid_t const,
                            fp_t * const, fp_t * const, fp_t * const);


/** @brief Calculate cloud optics using the parameterization described in
           https://doi.org/10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int hu_stamnes_1993(fp_t const liquid_water_path, /**< Liquid water path [g m-3].*/
                           fp_t const droplet_equivalent_radius, /**< Droplet equivalent radius [microns].*/
                           SpectralGrid_t const grid, /**< Spectral grid.*/
                           fp_t * const optical_depth, /**< Optical depth (wavenumber).*/
                           fp_t * const single_scatter_albedo, /**< Single-scatter albedo (wavenumber).*/
                           fp_t * const asymmetry_factor /**< Asymmetry factor (wavenumber).*/
                          );


/** @brief Calculate cloud optics using the parameterization described in
           https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int slingo_1989(fp_t const liquid_water_path, /**< Liquid water path [g m-3].*/
                       fp_t const droplet_equivalent_radius, /**< Droplet equivalent radius [microns].*/
                       SpectralGrid_t const grid, /**< Spectral grid.*/
                       fp_t * const optical_depth, /**< Optical depth (wavenumber).*/
                       fp_t * const single_scatter_albedo, /**< Single-scatter albedo (wavenumber).*/
                       fp_t * const asymmetry_factor /**< Asymmetry factor (wavenumber).*/
                      );


/** @brief Calculate liquid cloud optics.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int cloud_optics(Optics_t * const optics, /**< Optics.*/
                        fp_t * const liquid_water_path, /**< Liquid water path [g m-3] (layer).*/
                        fp_t * const droplet_equivalent_radius, /**< Droplet equivalent radius [microns] (layer).*/
                        LiquidCloud_t parameterization /**< Parameterization.*/
                       );


#endif
