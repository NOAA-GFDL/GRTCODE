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

#ifndef SET_PARSEHITRANFILE_H_
#define SET_PARSEHITRANFILE_H_

#include <stdint.h>

#include "myreal.h"  /* defines REAL_t */

#ifdef __NVCC__
#include "cudaHelpers.cuh"
#endif

/* BEGIN HITRAN FORMATTING FORWARD CONFIG */
typedef union HITRAN2012_vals_t {
  void* nil;
  uint8_t i;
  float f;
  double d;
} HITRAN2012_vals_t;

typedef enum HITRAN2012_cols_t
{
  mol_c,
  iso_c,
  Vnn_c,
  Snn_c,
  A_c,
  Yair_c,
  Yself_c,
  Elo_c,
  n_c,
  del_c,
  Vu_c,
  Vl_c,
  Qu_c,
  Ql_c,
  Ierr_c,
  Iref_c,
  flag_c,
  gu_c,
  gl_c,
  NCOLS
} HITRAN2012_col_t;

typedef enum LookupCast_t {
  NIL,
  UI8,
  F32,
  F64
} LookupCast_t;

typedef enum RefLinePtrIdx_t {
  mol_pidx,
  iso_pidx,
  Vnn_pidx,
  Snn_ref_pidx,
  Yair_pidx,
  Yself_pidx,
  En_pidx,
  n_pidx,
  d_pidx
} RefLinePtrIdx_t;

#ifndef __NVCC__
const unsigned int HITRAN2012_fmt[NCOLS][2];
const unsigned int HITRAN2012_recordLen; /* 160 plus  */
const unsigned int HITRAN2012_pad; /* trailing newline char,plus sentinel */
#endif

/* END HITRAN FORMATTING FORWARD CONFIG */

typedef struct SingleRefLine_t {
  uint8_t iso;
  REAL_t Vnn;
  REAL_t Snn_ref;
  float Yair;
  float Yself;
  float En;
  float n;
  float d;
} SingleRefLine_t ;

typedef struct RefLinePtrs_t {
  uint8_t* iso;
  REAL_t* Vnn;
  REAL_t* Snn_ref;
  float* Yair;
  float* Yself;
  float* En;
  float* n;
  float* d;
  uint8_t mol;
  unsigned int nLines;
}  RefLinePtrs_t;

typedef struct RefLine_flags_t {
  /* unsigned int memset_const: 1; */
  /* unsigned int memset_const_addr: 1; */
  /* unsigned int f_is_writable : 1 ; */
  unsigned int cumemset_host_flags;
  unsigned int host : 1 ;
  unsigned int device : 1;
} RefLine_flags_t;

/* exposed prototypes */
RefLinePtrs_t parseHITRANfile( char* , RefLine_flags_t, REAL_t loWn, REAL_t hiWn);
RefLinePtrs_t allocHost( unsigned int, RefLine_flags_t) ;
void reallocRefLines( RefLinePtrs_t* , RefLine_flags_t , RefLine_flags_t);
#ifdef __NVCC__
RefLinePtrs_t allocDevice( unsigned int );
void freeDevice( RefLinePtrs_t );
#endif
void freeHost( RefLinePtrs_t , RefLine_flags_t);


#endif
