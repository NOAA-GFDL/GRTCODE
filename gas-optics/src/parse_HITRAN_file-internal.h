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

#ifndef PARSE_HITRAN_FILE_INTERNAL_H_
#define PARSE_HITRAN_FILE_INTERNAL_H_

#include "grtcode_utilities.h"
#include "parse_HITRAN_file.h"


/** @brief Free memory used for molecular line parameters.
    @return GRTCODE_SUCCESS or an error code.*/
int free_line_params(LineParams_t * const line_params /**< Molecular line parameters.*/
                    );


/** @brief Read molecular line parameters from a HITRAN file.
    @return GRTCODE_SUCCESS or an error code.*/
int parse_hitran_file(LineParams_t * const line_params, /**< Molecular line paraemeters.*/
                      char const * const filepath, /**< Path to the HITRAN file.*/
                      int const mol_id, /**< Molecule id.*/
                      double const w0, /**< Lower bound [cm-1] of the spectral range being considered.*/
                      double const wn, /**< Upper bound [cm-1] of the spectral range being considered.*/
                      Device_t const device /**< Device id.*/
                     );


#endif
