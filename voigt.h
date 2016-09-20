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

#ifndef _SET_VOIGT_H_
#define _SET_VOIGT_H_

#include "myreal.h"

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t eta(REAL_t lorFWHM, REAL_t gauFWHM);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauFWHM(REAL_t T, REAL_t M, REAL_t v0);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauAlphad(REAL_t T, REAL_t M, REAL_t v0);

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauKernel( REAL_t v, REAL_t v0, REAL_t alphad );
      
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t pseudoVoigt(REAL_t eta, REAL_t lory, REAL_t gauy);


#endif
