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

#include <math.h>

#include "voigt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif

/* implimentation of method: 
   Ida, T and Ando, M and Toraya, H (2000)
   "Extended pseudo-Voigt function for approximating the Voigt profile".
   Journal of Applied Crystallography 33 (6): 1311-1316.
   doi:10.1107/s0021889800010219

   Authors claims accurate to 1%
*/


#ifdef __NVCC__
__host__ __device__
#endif
static REAL_t f(REAL_t lorFWHM, REAL_t gauFWHM){
  REAL_t f;
  /* gaussian powers */
  const REAL_t gauFWHM2 = gauFWHM*gauFWHM;
  const REAL_t gauFWHM3 = gauFWHM2*gauFWHM;
  const REAL_t gauFWHM4 = gauFWHM2*gauFWHM2;
  const REAL_t gauFWHM5 = gauFWHM3*gauFWHM2;
  /* lorentzian powers */
  const REAL_t lorFWHM2 = lorFWHM*lorFWHM;
  const REAL_t lorFWHM3 = lorFWHM2*lorFWHM;
  const REAL_t lorFWHM4 = lorFWHM2*lorFWHM2;
  const REAL_t lorFWHM5 = lorFWHM3*lorFWHM2;
  /* coefs */
  const REAL_t c1 = 2.69269;
  const REAL_t c2 = 2.42843;
  const REAL_t c3 = 4.47163;
  const REAL_t c4 = 0.07842;

  f = gauFWHM5 + c1*gauFWHM4*lorFWHM + c2*gauFWHM3*lorFWHM2
      + c3*gauFWHM2*lorFWHM3 + c4*gauFWHM*lorFWHM4 + lorFWHM5;

  f = powf(f,0.2f);
  
  return f;
}
  
#ifdef __NVCC__
__host__ __device__
#endif
static REAL_t eta_(REAL_t lorFWHM, REAL_t f){
  REAL_t eta;
  /* reuse ratio */
  const REAL_t flf = lorFWHM/f;
  /* coefs */
  const REAL_t c1 = 1.36603;
  const REAL_t c2 = -0.47719;
  const REAL_t c3 = 0.11116;

  eta = flf*(c1 + c2*flf + c3*flf*flf);

  return eta;
}
  
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t eta(REAL_t lorFWHM, REAL_t gauFWHM){
  return eta_(lorFWHM, f(lorFWHM, gauFWHM) );
}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauAlphad(REAL_t T, REAL_t M, REAL_t v0){
  /* T  Kelvin
     M  Molar Mass
     v0 line center
  */
  /* const REAL_t R = 8.3144598f; /\* gas constant J mol-1 k-1 *\/ */
  /* const REAL_t R = 82.057338;  /\* gas constant  cm3 atm K-1 mol-1 *\/ */
  /* const REAL_t R = 0.082057338;  /\* gas constant  L atm K-1 mol-1 *\/ */
  /* const REAL_t alphad = v0 * sqrtf( (2*R*T)/M ); */
  /*  From 'Petty' */
  const REAL_t m = M/6.023E23;
  /* const REAL_t c =  299792458;  /\* m/s SI*\/ */
  /* const REAL_t kb = 1.38064852-E23  /\* SI  J/K = m2·kg/(s2·K)*\/ */
  /* CGS units! */
  const REAL_t c = 2.99792458E10;  /* cm/s */
  const REAL_t kb = 1.380658E-16;  /* erg/K-1 = g·cm2/s2 */
  const REAL_t alphad = v0*sqrtf( (2*kb*T)/(m*c*c ));

  return alphad;

}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauFWHM(REAL_t T, REAL_t M, REAL_t v0){
  /* T  Kelvin
     M  Molar Mass
     v0 line center
  */
  const REAL_t sqrtln2 = 0.83255461115f;
  const REAL_t HWHM = sqrtln2*gauAlphad(T, M, v0);
  return 2*HWHM;
}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t gauKernel( REAL_t v, REAL_t v0, REAL_t alphad ){
  const REAL_t x = (v-v0);
  return expf(-x*x/(alphad*alphad)) / (alphad*sqrtf(M_PI)) ;
}

      
#ifdef __NVCC__
__host__ __device__
#endif
REAL_t pseudoVoigt(REAL_t eta, REAL_t lory, REAL_t gauy){
  return eta*lory + (1.-eta)*gauy;
}
