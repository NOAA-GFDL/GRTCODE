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

#ifndef TIPS2017_H_
#define TIPS2017_H_

#include "debug.h"
#include "grtcode_utilities.h"


/** @brief Load look-up tables onto the GPU.
    @return GRTCODE_SUCCESS or an error code.*/
int inittips_d(void);


/** @brief Calculate the total partition function.
    @return Total partition function.*/
HOST DEVICE fp_t Q(int const mol_id, /**< Molecule id.*/
                   fp_t const T, /**< Temperature [K].*/
                   int const iso /**< Isotopologue id.*/
                  );


#endif
