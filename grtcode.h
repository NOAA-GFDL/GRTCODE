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

#ifndef SET_ABSORPTION_COEF_H_
#define SET_ABSORPTION_COEF_H_

#include <stdint.h>
#include "parseNetcdfRadiation.h"

#include "myreal.h"

/* opticsbufptrs_t */
typedef struct OpticsBufPtrs_t {
  REAL_t* Gam;
  REAL_t* PShift;
  REAL_t* S;
} OpticsBufPtrs_t;



#ifdef __NVCC__
#ifdef __cplusplus
extern "C"
#endif
__host__ int device_launch(const unsigned int numMols,
                           char**  molFnames,
                           const REAL_t loWn,
                           const REAL_t hiWn,
                           const unsigned int nF,
                           const REAL_t resolution,
                           const unsigned int wingBreadth,
                           radiationOutputFields_t* atmosData,
                           REAL_t* const out);
#endif

#ifdef __NVCC__
#ifdef __cplusplus
extern "C"
#endif
int host_launch(const unsigned int numMols,
                char**  molFnames,
                const REAL_t loWn,
                const REAL_t hiWn,
                const unsigned int nF,
                const REAL_t resolution,
                const unsigned int wingBreadth,
                radiationOutputFields_t* atmosData,
                REAL_t* const out);
#endif



#endif
