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

#ifndef MOLECULES_INTERNAL_H_
#define MOLECULES_INTERNAL_H_

#include "grtcode_utilities.h"
#include "molecules.h"


/** @brief Initialize a molecule.  Read in its line parameters from the
           input HITRAN database file.
    @return GRTCODE_SUCCESS or an error code.*/
int create_molecule(Molecule_t * const mol, /**< Molecule.*/
                    int const id, /**< HITRAN molecule id.*/
                    char const * const hitran_path, /**< Path to HITRAN database file.*/
                    double const w0, /**< Spectral lower bound [cm-1].*/
                    double const wn, /**< Spectral upper bound [cm-1].*/
                    int const num_layers, /**< Number of atmospheric layers.*/
                    Device_t const device /**< Device id.*/
                   );


/** @brief Free memory stored by the molecule object.
    @return GRTCODE_SUCCESS or an error code.*/
int free_molecule(Molecule_t * const mol /**< Molecule.*/
                 );


/** @brief Convert from a molecule id to an array index.
    @return GRTCODE_SUCCESS or an error code.*/
int molecule_hash(int const mol_id, /**< Molecule id.*/
                  int * const hash /**< Array index.*/
                 );


#endif
