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

#ifndef SET_PARSENETCDFRADIATION_H_
#define SET_PARSENETCDFRADIATION_H_

#include "myreal.h"

typedef struct radiationInputFields_t {
  float* RH2O;
  float* QO3;
  float* DPFLUX;
  float* PRESSM;
  float* TEMP;
  float* DELTAZ;
  size_t nlat;
  size_t nlon;
  size_t npfull;
  size_t nphalf;
  size_t ntime;
} radiationInputFields_t;

typedef struct radiationOutputFields_t {
  REAL_t* N;
  REAL_t* P;
  REAL_t* T;
  REAL_t* DELTAZ;
  REAL_t* PS;
  size_t nlat;
  size_t nlon;
  size_t npfull;
  size_t nphalf;
  size_t ntime;
} radiationOutputFields_t;

int radiationInputFieldsMalloc(radiationInputFields_t* in);
int radiationOutputFieldsMalloc(radiationOutputFields_t* out);
int radiationInputFieldsFree(radiationInputFields_t* in);
int radiationOutputFieldsFree(radiationOutputFields_t* out);
int readInputFieldsFromFile(char fname[], radiationInputFields_t* in);
REAL_t getNumberDensity( const REAL_t rh2o , const REAL_t dpflux, const int hitranMolId );
REAL_t getPartialPres(REAL_t rh2o, REAL_t pressm, const REAL_t molarMass);
int setOutputFields(radiationInputFields_t *in, radiationOutputFields_t *out);
int getAndSetAtmosFieldsFromFile(char fname[], radiationOutputFields_t* out);

int test(char fname[]);

#endif
