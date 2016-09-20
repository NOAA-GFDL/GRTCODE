/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef SET_CONTINUUM_H_
#define SET_CONTINUUM_H_

#include "myreal.h"

void get_CTM(REAL_t* const OPT_CTM,
             REAL_t const * const T,
             REAL_t const * const P,
             REAL_t const * const HGT,
             REAL_t const * const H2O,
             const size_t NLINES,
             const size_t NLAYERS);

#endif
