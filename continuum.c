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
#include <math.h>
#include "continuum.h"
#include "myreal.h"

static const int MAXCHARSPERLINE=128;

void parseCKD(const char fname[],
              REAL_t** AryPtr,
              const int maxwavenum) /* zero based array, logic contingent on CTM arith being 1 wavenum res */
{  
  FILE* F;
  int res;
  char line[MAXCHARSPERLINE];
  double v0;
  double v1;
        
  /* if not malloc'd */
  if(*AryPtr==NULL)
  {
    /* malloc zeros */
    *AryPtr = (REAL_t*)calloc(maxwavenum,sizeof(**AryPtr));
    if(*AryPtr==NULL)
    {
      /* check malloc */
      fprintf(stderr,"Malloc AryPtr Failed.  Aborting\n");
      exit(EXIT_SUCCESS);
    }
  }

  /* open */
  F = fopen(fname, "r");
  if( F==NULL )
  {
    /* check open */
    fprintf(stderr,"Open %s Failed.  Aborting\n", fname);
    exit(EXIT_SUCCESS);
  }
  
  while( fgets(line, MAXCHARSPERLINE, F) !=NULL)
  {
    sscanf(line, "%lf %lf", &v0, &v1);
      /* 2nd col ignored */  /* 3rd col ignored */
    res = ((int)v0) -1 ;
    
    /* printf("%f\t%f \t\t~~~>\t\t",v0,v1); */
    /* printf("%d\t%g \n",res,((REAL_t)v1)); */

    /* writes second col value at index prescribed by first col */
    /* assumes Wavenumbers start at 1, and array is zero based index */
    if(res<maxwavenum)  /* skip this wavenumber */
    {
      (*AryPtr)[ res ] = ((REAL_t)v1);
    }
  }

  res = fclose(F);
  if(res!=0)
  {
    fprintf(stderr, "Attempt to close file %s previously opened with handle at %p failed with %d\n\tAborting.\n", fname, F, res);
  }

}

/* void ctm_free(REAL_t** CS, */
/*               REAL_t** CF, */
/*               REAL_t** T0, */
/*               REAL_t** T0F)     */
/* { */
/*   free(*CS); */
/*   *CS = NULL; */
/*   free(*CF); */
/*   *CF = NULL; */
/*   free(*T0); */
/*   *T0 = NULL; */
/*   free(*T0F); */
/*   *T0F = NULL;   */
/* } */

void ctm_read(REAL_t** CS,
              REAL_t** CF,
              REAL_t** T0,
              REAL_t** T0F,
              int maxwavenum)    
{
  const char testCKDF_f[] = "continuum/CKDF.ppp";            /* T0F */
  const char testCKDS_f[] = "continuum/CKDS.ppp";            /* T0  */
  const char testMTCKDF_f[] = "continuum/296MTCKD25_F.ccf";  /* CF  */
  const char testMTCKDS_f[] = "continuum/296MTCKD25_S.ccf";  /* CS  */

  parseCKD(testCKDF_f, T0F, maxwavenum);
  parseCKD(testCKDS_f, T0, maxwavenum);
  parseCKD(testMTCKDF_f, CF, maxwavenum);
  parseCKD(testMTCKDS_f, CS, maxwavenum);
  
  return;
}

void get_PD(REAL_t* const PD, REAL_t const * const P, REAL_t const * const H2O, const size_t nlyr)
{
  unsigned int n;
  for(n = 0; n < nlyr ; ++n)
  {
    /* PD[n] = ((REAL_t)100.) * H2O[n] * P[n]/((REAL_t)1E6); original was written for P in mb */
    PD[n] = ((REAL_t)100.) * H2O[n] * (1013.25*P[n])/((REAL_t)1E6);
  }
}


void get_PF(REAL_t* const  PF, REAL_t const * const P, REAL_t const * const PD, const size_t nlyr)
{
  unsigned int n;
  for(n = 0; n < nlyr ; ++n)
  {
    /* PF[n] = ((REAL_t)100.) * P[n] - PD[n];  original was written for P in mb */
    PF[n] = ((REAL_t)100.) * (1013.25*P[n]) - PD[n];
  }
}


void get_TF_TD(REAL_t* const TF,  /* out */
               REAL_t* const TD,  /* out */
               REAL_t const * const PD,  /* from get_PD, requres P, H2O, N */
               REAL_t const * const PF,  /* from get_PF, requires P, PD, N */
               REAL_t const * const CS,  /* ctm read */
               REAL_t const * const CF,  /* ctm read */
               REAL_t const * const HGT, /* model layer height*/
               REAL_t const * const T,   /* model temp */
               REAL_t const * const T0,  /* ctm read */
               REAL_t const * const T0F, /* ctm read */
               const size_t nlines,
               const size_t nlyr)
{  
  unsigned int lyr;
  unsigned int lp;
  const REAL_t c1 = 1E3;
  const REAL_t c2 = 101325.;
  const REAL_t c3 = 1.3806E-19;
  const REAL_t ct = 296.0;
  
  for(lyr = 0; lyr < nlyr; ++lyr)
  {
    for(lp = 0; lp < nlines; ++lp)
    {      
      TF[lyr*nlines + lp ] = ( PD[lyr] * PF[lyr] * CF[lp] * HGT[lyr] * c1 /
                  (c2 * c3 * T[lyr] ) ) *
          ( ct / T[lyr] ) * exp( T0F[lp] * (ct - T[lyr]) ) ;
      
      TD[lyr*nlines + lp ] = ( (PD[lyr]*PD[lyr]) * CS[lp] * HGT[lyr] * c1 /
                  (c2 * c3 * T[lyr] ) ) *
          ( ct / T[lyr] ) * exp( T0[lp] * (ct - T[lyr]) ) ;              
    }  /* for line */
  }    /* for layer */
}
            

void get_CTM(REAL_t* const OPT_CTM,
             REAL_t const * const T,
             REAL_t const * const P,
             REAL_t const * const HGT,
             REAL_t const * const H2O,
             const size_t NLINES,
             const size_t NLAYERS)
{
  unsigned int n;
  REAL_t* CS = NULL;
  REAL_t* CF = NULL;
  REAL_t* T0 = NULL;
  REAL_t* T0F = NULL;
  REAL_t* IBUF = NULL;
  REAL_t* PD = NULL;
  REAL_t* PF = NULL;
  REAL_t* PBUF = NULL;
  REAL_t* TD = NULL;
  REAL_t* TF = NULL;
  REAL_t* TBUF = NULL;

  /* later for gpu improvement,
   * we almost certainly want this to be array of float4,
   * not 4 arrays of floats... */
  
  IBUF = (REAL_t*)calloc(4*NLINES,sizeof(REAL_t));
  if(IBUF == NULL)
  {
    fprintf(stderr,"IBUF malloc failed");
    exit(EXIT_FAILURE);
  }
  CS = &(IBUF[0]);
  CF = &(IBUF[1*NLINES]);
  T0 = &(IBUF[2*NLINES]);
  T0F = &(IBUF[3*NLINES]);

  PBUF = (REAL_t*)calloc(2*NLAYERS,sizeof(REAL_t));
  if(PBUF == NULL)
  {
    fprintf(stderr,"PBUF malloc failed");
    exit(EXIT_FAILURE);
  }
  PD = &(PBUF[0]);
  PF = &(PBUF[NLAYERS]);

  TBUF = (REAL_t*)calloc(2*NLINES*NLAYERS,sizeof(REAL_t));
  if(TBUF == NULL)
  {
    fprintf(stderr,"TBUF malloc failed");
    exit(EXIT_FAILURE);
  }
  TF = &(TBUF[0]);
  TD = &(TBUF[NLINES*NLAYERS]);
  
  ctm_read(&CS, &CF, &T0, &T0F, NLINES);

  /* printf("%g %g %g\n", H2O[0],H2O[1],H2O[2]); */
  get_PD( PD, P, H2O, NLAYERS);
  get_PF( PF, P, PD, NLAYERS);
  /* for( n = 0; n < NLAYERS ; ++n ) */
  /* { */
  /*   printf("%g %g\n", PD[n], PF[n]);   */
  /* } */

  
  get_TF_TD( TF,  /* out */
             TD,  /* out */
             PD,  /* from get_PD, requres P, H2O, NLAYER */
             PF,  /* from get_PF, requires P, PD, NLAYER */
             CS,  /* ctm read */
             CF,  /* ctm read */
             HGT, /* model layer height*/
             T,   /* model temp */
             T0,  /* ctm read */
             T0F, /* ctm read */
             NLINES,
             NLAYERS);

  for( n = 0; n < NLINES*NLAYERS ; ++n )
  {
    /* sum self and foreign */
    OPT_CTM[n] += TF[n] + TD[n];
    /* /\* debug print *\/ */
    /* printf("DBG CTM: %g \n",TF[n] + TD[n]); */
  }
  
  free(IBUF);
  free(TBUF);
  free(PBUF);
  
}
  

#ifndef SKIPMAIN
/* tests */
int main(int argc, char** argv)
{
  int l;
  unsigned int n;
  const size_t N = 50000;
  const size_t NLYR = 3;
  
  REAL_t T[NLYR];
  REAL_t P[NLYR];
  REAL_t HGT[NLYR];
  REAL_t H2O[NLYR];
  REAL_t OPT[N*NLYR];
  for(n=0; n<NLYR; ++n)
  {
    T[n] = 276.;
    P[n] = 1.;
    HGT[n] = 1.;
    H2O[n] = 1.876E+04;
  }
  /* OPT[N*NLYR] = {0...}; */
  
  get_CTM(OPT,
          T,
          P,
          HGT,
          H2O,
          N,
          NLYR);

  for(n=0; n<N; ++n){
    for( l=0; l<NLYR; ++l){
      printf("%f ",OPT[l*N +n]);  /*  print layers as col */
    }
    printf("\n");  /* new line for each row */
  }
       
  exit(EXIT_SUCCESS);
}
#endif
