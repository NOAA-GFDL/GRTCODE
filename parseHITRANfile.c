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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "parseHITRANfile.h"

/* BEGIN HITRAN FORMATTING CONFIG */
const unsigned int HITRAN2012_fmt[NCOLS][2] = {
  {2 , UI8},
  {1 , UI8},
  {12, F64},
  {10, F64},
  {10, NIL},
  {5 , F32},
  {5,  F32},
  {10, F32},
  {4,  F32},
  {8,  F32},
  {15, NIL},
  {15, NIL},
  {15, NIL},
  {15, NIL},
  {6,  NIL},
  {12, NIL},
  {1,  NIL},
  {7,  NIL},
  {7,  NIL}
};
const unsigned int HITRAN2012_recordLen = 160;  /* 160 plus  */
const unsigned int HITRAN2012_pad = 2;  /* trailing newline char,plus sentinel */
/* END HITRAN FORMATTING CONFIG */


RefLinePtrs_t allocHost(unsigned int nLines, RefLine_flags_t flags) {
  const unsigned int cuFlags = flags.cumemset_host_flags;
  RefLinePtrs_t self;
  self.nLines = nLines;
  self.mol = 0;

  if (cuFlags != ((unsigned int)-1) ) {
#ifdef __NVCC__
    /* cudaHostAlloc */
    HANDLE_ERROR( cudaHostAlloc( &(self.iso) , nLines*sizeof(*(self.iso)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.Vnn) , nLines*sizeof(*(self.Vnn)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.Snn_ref) , nLines*sizeof(*(self.Snn_ref)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.Yair) , nLines*sizeof(*(self.Yair)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.Yself) , nLines*sizeof(*(self.Yself)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.En) , nLines*sizeof(*(self.En)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.n) , nLines*sizeof(*(self.n)), cuFlags) );
    HANDLE_ERROR( cudaHostAlloc( &(self.d) , nLines*sizeof(*(self.d)), cuFlags) );
#else
    fprintf(stderr,"You have not built for CUDA. This will end badly. Aborting.\n.");
    exit(1);
#endif
  }
  else
  {
    /* malloc */
    self.iso = ( uint8_t* )malloc( nLines*sizeof(*(self.iso)) ) ;
    self.Vnn = ( REAL_t* )malloc( nLines*sizeof(*(self.Vnn)) );
    self.Snn_ref = ( REAL_t* )malloc( nLines*sizeof(*(self.Snn_ref)) );
    self.Yair = ( float* )malloc( nLines*sizeof(*(self.Yair)) );
    self.Yself = ( float*)malloc( nLines*sizeof(*(self.Yself)) );
    self.En = ( float* )malloc( nLines*sizeof(*(self.En)) );
    self.n = ( float* )malloc( nLines*sizeof(*(self.n)) );
    self.d = ( float* )malloc( nLines*sizeof(*(self.d)) );

    if ( self.iso == NULL ||
         self.Vnn == NULL ||
         self.Snn_ref == NULL ||
         self.Yair == NULL ||
         self.Yself == NULL ||
         self.En == NULL ||
         self.n == NULL ||
         self.d == NULL	){
      fprintf(stderr,"allocHost failed to malloc\n");
      exit(1);
    }
  }
  return self;
}

#ifdef __NVCC__
RefLinePtrs_t allocDevice( unsigned int nLines ) {
  RefLinePtrs_t self;
  self.nLines = nLines;
  self.mol =0;
  
  /* cudaDeviceAlloc */
  HANDLE_ERROR( cudaMalloc( &(self.iso) , nLines*sizeof( *(self.iso)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.Vnn) , nLines*sizeof( *(self.Vnn)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.Snn_ref) , nLines*sizeof( *(self.Snn_ref)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.Yair) , nLines*sizeof( *(self.Yair)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.Yself) , nLines*sizeof( *(self.Yself)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.En) , nLines*sizeof( *(self.En)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.n) , nLines*sizeof( *(self.n)) ) );
  HANDLE_ERROR( cudaMalloc( &(self.d) , nLines*sizeof( *(self.d)) ) );

  return self;
}

void freeDevice(RefLinePtrs_t self){
  HANDLE_ERROR( cudaFree( self.iso ) );
  HANDLE_ERROR( cudaFree( self.Vnn ) );
  HANDLE_ERROR( cudaFree( self.Snn_ref ) );
  HANDLE_ERROR( cudaFree( self.Yair ) );
  HANDLE_ERROR( cudaFree( self.Yself ) );
  HANDLE_ERROR( cudaFree( self.En ) );
  HANDLE_ERROR( cudaFree( self.n ) );
  HANDLE_ERROR( cudaFree( self.d ) );
  return;
}
#endif

void freeHost(RefLinePtrs_t self, RefLine_flags_t flags){
  const unsigned int cuFlags = flags.cumemset_host_flags;

  if ( cuFlags != ((unsigned int)-1) ) {
#ifdef __NVCC__
    /* cudaHostAlloc */
    HANDLE_ERROR( cudaFreeHost( self.iso ) );
    HANDLE_ERROR( cudaFreeHost( self.Vnn ) );
    HANDLE_ERROR( cudaFreeHost( self.Snn_ref ) );
    HANDLE_ERROR( cudaFreeHost( self.Yair ) );
    HANDLE_ERROR( cudaFreeHost( self.Yself ) );
    HANDLE_ERROR( cudaFreeHost( self.En ) );
    HANDLE_ERROR( cudaFreeHost( self.n ) );
    HANDLE_ERROR( cudaFreeHost( self.d ) );
#else
    fprintf(stderr,"You have not built for CUDA. This will end badly. Aborting.\n.");
    exit(1);
#endif

  }
  else {
    free( self.iso );
    free( self.Vnn );
    free( self.Snn_ref );
    free( self.Yair );
    free( self.Yself );
    free( self.En );
    free( self.n );
    free( self.d );
  }
  return;
}

void reallocRefLines( RefLinePtrs_t* Lines, RefLine_flags_t old_flags, RefLine_flags_t new_flags){
  /* to realloc as minimal size and/or different flags */
  unsigned int i;
  /* shallow copy */
  RefLinePtrs_t old = *Lines;

  /* alloc a smaller array */
  RefLinePtrs_t newLines = allocHost(old.nLines , new_flags);

  /* deep copy */
  newLines.mol = old.mol;
  newLines.nLines = old.nLines;
  for(i=0 ; i<old.nLines ; ++i){
    newLines.iso[i] = old.iso[i];
    newLines.Vnn[i] = old.Vnn[i];
    newLines.Snn_ref[i] = old.Snn_ref[i];
    newLines.Yair[i] = old.Yair[i];
    newLines.Yself[i] = old.Yself[i];
    newLines.En[i] = old.En[i];
    newLines.n[i] = old.n[i];
    newLines.d[i] = old.d[i];
  }

  /* free old */
  freeHost(old, old_flags);

  /* update pointer */
  *Lines = newLines;

  return;
}


HITRAN2012_vals_t HITRAN2012_cast(const int col, char* sval){
  const LookupCast_t typ = (LookupCast_t)HITRAN2012_fmt[col][1];
  HITRAN2012_vals_t val={0};
  char* endptr;
  errno=0;
  switch (typ) {
    case NIL :
      val.nil= NULL;
      break;
    case UI8 :
      val.i = (uint8_t)( strtol(sval, &endptr,10) );
      if (val.i==0 && errno !=0 ){
        /* error */
        if (errno == EINVAL){
          fprintf(stderr,"strtol failed: Invalid Value.\nGiven %s \nLeaving tail of %s\nAborting.\n", sval, endptr);
          exit(1);
        }
        else if (errno == ERANGE){
          fprintf(stderr,"strtol reported out of range: \nGiven %s leaving tail of %s\nStoring %d and Continuing\n", sval, endptr,val.i);
        }
        else {
          fprintf(stderr,"strtol failed: Uknown Error, errno: %d .\nGiven %s \nLeaving tail of %s\nAborting.\n", errno, sval, endptr);
          exit(1);
        }
      }
      break;
    case F64:
    case F32:
      /* falls through from F64 */
      val.d = (strtod(sval, &endptr));
      if (val.d==0 && errno !=0 ){
        /* error */
        if (errno == EINVAL){
          fprintf(stderr,"strtod failed: Invalid Value.\nGiven %s \nLeaving tail of %s\nAborting.\n", sval, endptr);
          exit(1);
        }
        else if (errno == ERANGE){
          fprintf(stderr,"strtod reported out of range: \nGiven %s leaving tail of %s\nStoring %f and Continuing\n", sval, endptr,val.d);
        }
        else {
          fprintf(stderr,"strtod failed: Uknown Error, errno: %d .\nGiven %s \nLeaving tail of %s\nAborting.\n", errno, sval, endptr);
          exit(1);
        }
      }
      /* for floats, cast */
      if (typ==F32){
        val.f = val.d;  /* cast to float */
      }
      break;

    default :
      fprintf(stderr,"\nHITRAN2012_cast failed on col %d, LookupCast_t %d, sval: %s \n",HITRAN2012_fmt[col][0],HITRAN2012_fmt[col][1],sval);
      exit(1);
  }

  return val;
}



/*/\* for portability, though I haven't seen a machine without getline in a while...*\/ */
/* int getsLineFromFile(char* line, int maxLineLen, FILE *fp){ */
/*   if ( fgets(line,maxLineLen+1,fp) == NULL){ */
/*     return -1; */
/*   } */
/*   else { */
/*     return strnlen(line,maxLineLen); */
/*   } */
/* } */


RefLinePtrs_t parseHITRANfile( char fname[], RefLine_flags_t flags, REAL_t loWn, REAL_t hiWn){
  unsigned int n=0;
  size_t l=0;
  ssize_t ll=0;
  const unsigned int MAXLINE=163;
  char* buf = (char*)calloc(MAXLINE,sizeof(char));

  RefLinePtrs_t RefLines;

  /* parsing temps */
  int col;
  size_t offset;
  size_t len;
  char tmp[16];
  unsigned int valIdx;  /* begin */
  unsigned int t;
  
  /* open file */
  FILE *fp = fopen(fname,"r");
  if (fp == NULL){
    fprintf(stderr,"error opening file %s \n",fname);
    exit(1);
  }

  /* line count */
  /* while( (l=getsLineFromFile(buf, MAXLINE, fp)) != -1) */
  while ( (ll=getline(&buf,&l,fp)) != -1)
  {
    ++n;
  }
  fprintf(stderr,"Found %s has %d total lines. Parsing into struct.\n",fname,n);



  /* malloc */
  RefLines = allocHost(n,flags);

  /* parse */
  rewind(fp);
  n=0;
  while( (ll=getline( &buf, &l, fp)) != -1){
    if ( (ll-HITRAN2012_recordLen)>HITRAN2012_pad || ll<HITRAN2012_recordLen){
      if (ll > (HITRAN2012_recordLen + HITRAN2012_pad) ){
        fprintf(stderr,"\nFound bad record at line %d ( %lu exceeds max %d chars) in %s\n",n,l,MAXLINE,fname);
      }
      else{
        fprintf(stderr,"\nFound bad record (too short) at line %d in %s, len %lu \n",n,fname,l);
      }
      exit(1);
    }

    valIdx=0;  /* begin */
    offset=0;
    /* read format, place into arrays */
    for (col=0; col<NCOLS ; ++col){
      len = HITRAN2012_fmt[col][0];
      t =HITRAN2012_fmt[col][1];
      strncpy(tmp , &(buf[offset]), len);
      tmp[len]='\0';  /* sentinel */
      offset+=len;
      if ( t != NIL){
        switch( valIdx ){
          case mol_pidx:
            t = HITRAN2012_cast( col , tmp).i;
            if (n==0){
              RefLines.mol = t;
            }
            else if( RefLines.mol != t ){
              fprintf(stderr,"Molecule changed from %d to %d at line %d. Aborting\n",RefLines.mol, t, n);
              exit(1);
            }
            break;
          case iso_pidx :
            RefLines.iso[n] = HITRAN2012_cast( col , tmp).i;
            break;
          case Vnn_pidx:
            RefLines.Vnn[n] = (REAL_t)HITRAN2012_cast( col , tmp).d;
            break;
          case   Snn_ref_pidx:
            RefLines.Snn_ref[n] = (REAL_t)HITRAN2012_cast( col , tmp).d;
            break;
          case   Yair_pidx:
            RefLines.Yair[n] = HITRAN2012_cast( col , tmp).f;
            break;
          case   Yself_pidx:
            RefLines.Yself[n] = HITRAN2012_cast( col , tmp).f;
            break;
          case   En_pidx:
            RefLines.En[n] = HITRAN2012_cast( col , tmp).f;
            break;
          case   n_pidx:
            RefLines.n[n] = HITRAN2012_cast( col , tmp).f;
            break;
          case d_pidx:
            RefLines.d[n] = HITRAN2012_cast( col , tmp).f;
            break;
          default:
            fprintf(stderr,"\n Unkown RefLinePtrIdx_t=%d, abort \n",valIdx);
            exit(1);
        }
        ++valIdx;
      }
    }
    if ( (loWn < 0) && (hiWn < 0) ) {
      /* always contribute to index */
      ++n;
    }
    else { 			/* filter */
      if( (RefLines.Vnn[n] >= loWn) && (RefLines.Vnn[n] <= hiWn) ){
	/*  contribute to index */
	++n;
      }
      /* else pretend it doesn't exist */
    }
  }
  /* set nLines */
  RefLines.nLines = n;
  
  if (fclose(fp)){
    fprintf(stderr,"error closing file %s \n",fname);
    exit(1);
  }

  if ( (loWn >=0) || (hiWn>=0) ){
    fprintf(stderr, "Filtering lines in range [%f,%f] resulted in %d lines.\n", loWn, hiWn, n);
    reallocRefLines(&RefLines,flags,flags);
  }
  
  return RefLines;
}

#ifndef SKIPMAIN
int main(int argc, char* argv[]){

  RefLine_flags_t flags= {0,1,0}; /* cumalloc default, host=True, device=false */
  RefLinePtrs_t L = parseHITRANfile(argv[1],flags);

  return EXIT_SUCCESS;
}
#endif
