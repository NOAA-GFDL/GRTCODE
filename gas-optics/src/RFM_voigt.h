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

#ifndef RFM_VOIGT_H_
#define RFM_VOIGT_H_

#include "debug.h"
#include "grtcode_utilities.h"
#include "line_shape.h"


/** @brief Calculate the Voigt line shape.
    @return GRTCODE_SUCCESS or an error code.*/
HOST DEVICE int rfm_voigt_line_shape(LineShapeInputs_t const vals, /**< Line shape parameters.*/
                                     fp_t * const K /**< Line shape values (wavenumber).*/
                                    );


#endif
