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

#ifndef WATER_VAPOR_CONTINUUM_INTERNAL_H_
#define WATER_VAPOR_CONTINUUM_INTERNAL_H_

#include "grtcode_utilities.h"
#include "water_vapor_continuum.h"


/** @brief Water vapor continuum coefficients.*/
enum water_vapor_coefs
{
    MTCKD25_F296 = 0,
    MTCKD25_S296,
    CKDF,
    CKDS,
    NUM_COEFS
};


/** @brief Read in water vapor continuum coefficients.
    @return GRTCODE_SUCCESS or an error code.*/
int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc, /**< Continum coefficents.*/
                                    char const * const h2o_ctm_dir, /**< Directory containing the input file.*/
                                    SpectralGrid_t const grid, /**< Spectral grid.*/
                                    Device_t const device /**< Device id.*/
                                   );


/** @brief Free memory used for the continuum coefficients.
    @return GRTCODE_SUCCESS or an error code.*/
int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc /**< Continuum coefficients.*/
                                    );


#endif
