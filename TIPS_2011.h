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

#ifndef SET_TIPS_2011_H_
#define SET_TIPS_2011_H_

#include <stdint.h>

#ifdef __NVCC__
__host__ __device__
#endif
void QT(const uint8_t molNum, /* HITRAN molecule ID number */
        const REAL_t T,       /* temperature in K */
        const uint8_t iso,    /* isotope code (HITRAN INDEX) */
        float* const gsi,    /* state independent nuclear degeneracyfactor */
        REAL_t* const Qt);    /* Total Internal Partition Function */

#ifdef __NVCC__
__host__ int initTIPS_d();
#endif

int initTIPS();

#endif
