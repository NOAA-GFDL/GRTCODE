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
#include <assert.h>
#include <netcdf.h>
#include "parseNetcdfRadiation.h"

static const int NUM_MOLS = 7;

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define NCERR(e) {fprintf(stderr, "Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}


int readInputFieldsFromFile(char fname[], radiationInputFields_t* in){

  /* Error handling. */
  int retval;
  /* NC Handles */
  int ncid;
  int varid;
  int dimid;
  /* We will learn about the data file and store results in these
     program variables. */
  int ndims_in;
  int nvars_in;
  int ngatts_in;
  int unlimdimid_in;


  /* open the file. */
  if ( (retval = nc_open(fname, NC_NOWRITE, &ncid) ) ){    
    NCERR(retval);
  }
  /* inquire about the file */
  if ( ( retval = nc_inq(ncid,
                         &ndims_in,
                         &nvars_in,
                         &ngatts_in,
                         &unlimdimid_in ) ) ){
    NCERR(retval);
  }

  /* If we needed to test things regarding the file,
     they'd go here */
  /*  */

  /* get cordinate system bounds from nc dimensions */
  /* lat */
  if ( ( retval = nc_inq_dimid(ncid, "lat", &dimid ) ) ){
    NCERR(retval);
  }
  if ((retval = nc_inq_dimlen(ncid, dimid, &in->nlat) ) ){
    NCERR(retval);
  }
  /* lon */
  if ( ( retval = nc_inq_dimid(ncid, "lon", &dimid ) ) ){
    NCERR(retval);
  }
  if ((retval = nc_inq_dimlen(ncid, dimid, &in->nlon) ) ){
    NCERR(retval);
  }
  /* pfull */
  if ( ( retval = nc_inq_dimid(ncid, "pfull", &dimid ) ) ){
    NCERR(retval);
  }
  if ((retval = nc_inq_dimlen(ncid, dimid, &in->npfull) ) ){
    NCERR(retval);
  }
  /* phalf */
  if ( ( retval = nc_inq_dimid(ncid, "phalf", &dimid ) ) ){
    NCERR(retval);
  }
  if ((retval = nc_inq_dimlen(ncid, dimid, &in->nphalf) ) ){
    NCERR(retval);
  }
  /* time */
  if ( ( retval = nc_inq_dimid(ncid, "time", &dimid ) ) ){
    NCERR(retval);
  }
  if ((retval = nc_inq_dimlen(ncid, dimid, &in->ntime) ) ){
    NCERR(retval);
  }  
  /* sanity check cell vs skeleton */
  assert( in->nphalf == in->npfull+1);  /* should have one more level than layer... */
  /* end cords */

  /* make some room */
  radiationInputFieldsMalloc(in);
        
  
  /* get the varids for the variable we care about */
  /* read the data for the variable we care  about */
  /* RH2O */
  if ( ( retval = nc_inq_varid(ncid, "rh2o", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->RH2O ) ) ){
    NCERR(retval);
  }
  /* qo3 */
  if ( ( retval = nc_inq_varid(ncid, "qo3", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->QO3 ) ) ){
    NCERR(retval);
  }
  /* DPFLUX */
  if ( ( retval = nc_inq_varid(ncid, "dpflux", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->DPFLUX ) ) ){
    NCERR(retval);
  }
  /* PRESSM */
  if ( ( retval = nc_inq_varid(ncid, "pressm", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->PRESSM ) ) ){
    NCERR(retval);
  }
  /* TEMP */
  if ( ( retval = nc_inq_varid(ncid, "temp", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->TEMP ) ) ){
    NCERR(retval);
  }
  /* DELTAZ */
  if ( ( retval = nc_inq_varid(ncid, "deltaz", &varid) ) ){
    NCERR(retval);
  }
  if ( ( retval = nc_get_var_float(ncid, varid, in->DELTAZ ) ) ){
    NCERR(retval);
  }

  
  return EXIT_SUCCESS;
}


int radiationInputFieldsMalloc(radiationInputFields_t* in){
  /* ASSUMES DIMS ALREADY PACKED IN STRUCT */
  /* n elements in a pfull field */
  size_t pfulllen = in->ntime * in->npfull * in->nlat * in->nlon;
  /* n elements in a phalf field */
  size_t phalflen = in->ntime * in->nphalf * in->nlat * in->nlon;
  /* mallocs */
  in->RH2O   = (float*)malloc( pfulllen * sizeof(*(in->RH2O)) );
  in->QO3    = (float*)malloc( pfulllen * sizeof(*(in->QO3)) );
  in->DPFLUX = (float*)malloc( pfulllen * sizeof(*(in->DPFLUX)) );
  in->PRESSM = (float*)malloc( pfulllen * sizeof(*(in->PRESSM)) );
  in->TEMP   = (float*)malloc( pfulllen * sizeof(*(in->TEMP)) );
  in->DELTAZ = (float*)malloc( phalflen * sizeof(*(in->DELTAZ)) );
  if ( in->RH2O==NULL || in->QO3==NULL || in->DPFLUX==NULL || in->PRESSM==NULL || in->TEMP==NULL || in->DELTAZ==NULL)
  {
    fprintf(stderr, "radiationInputFieldsMalloc FAILED!\n");
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}

int radiationInputFieldsFree(radiationInputFields_t* in){
  /* really should init these null, check isnull when malloc, and check notnull when free,
     but we're using netcdf anyway, so I digress*/
  free(in->RH2O);
  free(in->QO3);
  free(in->DPFLUX);
  free(in->PRESSM);
  free(in->TEMP);
  free(in->DELTAZ);
  return EXIT_SUCCESS;
}  

int radiationOutputFieldsFree(radiationOutputFields_t* out){
  /* really should init these null, check isnull when malloc, and check notnull when free,
     but we're using netcdf anyway, so I digress*/
  free(out->N);
  free(out->P);
  free(out->T);
  free(out->DELTAZ);
  free(out->PS);
  return EXIT_SUCCESS;
}  


int radiationOutputFieldsMalloc(radiationOutputFields_t* out){
  /* ASSUMES DIMS ALREADY PACKED IN STRUCT */
  /* n elements in a pfull field */
  size_t pfulllen = out->ntime * out->npfull * out->nlat * out->nlon;
  /* n elements in a phalf field */
  /* size_t phalflen = out->ntime * out->nphalf * out->nlat * out->nlon; */
  /* mallocs */
  out->N = (REAL_t*)malloc( NUM_MOLS * pfulllen * sizeof(*(out->N)) );
  out->P = (REAL_t*)malloc( pfulllen * sizeof(*(out->P)) );
  out->T = (REAL_t*)malloc( pfulllen * sizeof(*(out->T)) );
  out->DELTAZ = (REAL_t*)malloc( pfulllen * sizeof(*(out->DELTAZ)) );
  /* currently supports first seven hitran molecules */
  out->PS = (REAL_t*)malloc( NUM_MOLS * pfulllen * sizeof(*(out->PS)) );
  if ( out->N==NULL || out->P==NULL || out->T==NULL || out->DELTAZ==NULL || out->PS==NULL){
    fprintf(stderr, "radiationOutputFieldsMalloc FAILED!\n");
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}

/* REAL_t getNumberDensity( const REAL_t ratio, const REAL_t dpflux , const int hitranMolId){ */
/*   REAL_t coef; */
/*   const REAL_t g = 9.81; */
/*   switch(hitranMolId) */
/*   { */
/*     case 1:           /\* H2O *\/ */
/*       coef = 3.34E21;  /\* per David's Email *\/ */
/*       break; */
/*     /\* case 2:  /\\* CO2 *\\/ *\/ */
/*     /\*   coef = 3.34E21; *\/ */
/*     /\*   break; *\/ */
/*     case 3: */
/*       coef = 3.34E21; */
/*       break; */
/*     /\* case 4: *\/ */
/*     /\*   coef = 3.34E21; *\/ */
/*     /\*   break; *\/ */
/*     /\* case 5: *\/ */
/*     /\*   coef = 3.34E21; *\/ */
/*     /\*   break; *\/ */
/*     /\* case 6: *\/ */
/*     /\*   coef = 3.34E21; *\/ */
/*     /\*   break; *\/ */
/*     /\* case 7: *\/ */
/*     /\*   coef = 3.34E21; *\/ */
/*     /\*   break; *\/ */
/*     default: */
/*       fprintf(stderr,"Conversion of HitranMolId %d is not currently implimented\n",hitranMolId); */
/*       exit(EXIT_FAILURE); */
/*   }   */
/*   return ratio * dpflux * coef / g ; */
/* } */


REAL_t getPartialPres(const REAL_t ratio, const REAL_t pressm, const REAL_t molarMass)
{
  const REAL_t molarMassDryAir= 29;
  
  return ratio * ( molarMassDryAir / molarMass ) * pressm ;
  
}


int setOutputFields(radiationInputFields_t *in, radiationOutputFields_t *out){
  /* Indexing */
  size_t ifoffset;  /* input data is time,pfull,lat,lon */
  size_t ofoffset; /* would prefer output data as time,lat,lon,pfull */
  size_t ihoffset;  /* input data is time,phalf(pfull+1),lat,lon */
  size_t h2ooffset;
  size_t o3offset;
  
  /* input data is time,pfull,lat,lon */
  size_t i;
  size_t j;
  size_t k;
  size_t t;

  /* copy meta data from Input Fields */ 
  out->nlat   = in->nlat;
  out->nlon   = in->nlon;
  out->npfull = in->npfull;
  out->nphalf = in->nphalf;
  out->ntime  = in->ntime;
  /* make some room */
  radiationOutputFieldsMalloc(out);

  /* reshape, I thought netcdf did this stuff, meh */
  for(t=0; t<in->ntime; ++t){
    for(k=0; k<in->npfull; ++k){
      for(i=0; i<in->nlat; ++i){
        for(j=0; j<in->nlon; ++j){
          /* found input and output offsets */
          ifoffset=t*(in->npfull*in->nlat*in->nlon) + k*(in->nlat*in->nlon) + i*(in->nlon) +j;
          ofoffset=t*(in->nlat*in->nlon*in->npfull) + i*(in->nlon*in->npfull) + j*(in->npfull) +k;
          ihoffset=t*(in->nphalf*in->nlat*in->nlon) + k*(in->nlat*in->nlon) + i*(in->nlon) +j;
          h2ooffset = t*(in->nlat*in->nlon*in->npfull*NUM_MOLS)
              + i*(in->nlon*in->npfull*NUM_MOLS)
              + j*(in->npfull*NUM_MOLS)
              + 0*in->npfull + k ;
          o3offset = t*(in->nlat*in->nlon*in->npfull*NUM_MOLS)
              + i*(in->nlon*in->npfull*NUM_MOLS)
              + j*(in->npfull*NUM_MOLS)
              + 2*in->npfull + k;
          /* pack the data into out */
          /* unit conv */
          out->P[ofoffset] = ((REAL_t)in->PRESSM[ifoffset])*9.86923E-6;  /* converts Pa to atm */
          out->T[ofoffset] = ((REAL_t)in->TEMP[ifoffset]);               /* Kelvin */
          /* calcs */
          /* if: method is already packing level path integral (ie, dpflux, cough) then: */
          /* out->DELTAZ[ofoffset] = ((REAL_t)1); */
          /* else: */
          out->DELTAZ[ofoffset] = in->DELTAZ[ihoffset] * 100.;  /* m to cm */
          
          out->PS[h2ooffset] = getPartialPres( ((REAL_t)in->RH2O[ifoffset]) , ((REAL_t)in->PRESSM[ifoffset]), 18. ) * 9.86923E-6 ; /* Ps, correction to atm */
          /* o3 is third hitran mol */
          out->PS[o3offset] = getPartialPres( ((REAL_t)in->QO3[ifoffset]) , ((REAL_t)in->PRESSM[ifoffset]), 48. ) * 9.86923E-6 ; /* Ps, correction to atm */
          
        }                                       /* t */
      }                                         /* k */
    }                                           /* lat */
  }                                             /* lon */

  /* don't think we need these anymore */
  radiationInputFieldsFree(in);
  
  return EXIT_SUCCESS;

}

int getAndSetAtmosFieldsFromFile(char fname[], radiationOutputFields_t* out){
  radiationInputFields_t in;

  if(fname==NULL || strlen(fname)==0 )
  {
    fprintf(stderr, "Invalid atmos file, did you specify an atmos file?!\n");
    exit(EXIT_FAILURE);
  }
  
  printf("Attempting to open and input from %s. \n",fname);
  readInputFieldsFromFile(fname, &in);

  printf("Attempting to reshape, compute, and set output.\n");
  setOutputFields(&in, out);

  return EXIT_SUCCESS;
}
  

int test(char fname[]){
  radiationOutputFields_t out;
  unsigned int i; /* = 0; */
  unsigned int j; /* = 0; */
  unsigned int k; /* = 47; */
  unsigned int t;

  unsigned int offset;

  getAndSetAtmosFieldsFromFile(fname, &out);

  printf("\nTest Results:\n\n");
  for(t=0; t<out.ntime; ++t){
    for(i=0; i<out.nlat; ++i){
      for(j=0; j<out.nlon; ++j){
        for(k=0; k<out.npfull; ++k){
          offset = t*(out.nlat*out.nlon*out.npfull) + i*(out.nlon*out.npfull) + j*(out.npfull) +k;
          printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.N[offset]);
        }
      }
    }
  }

  radiationOutputFieldsFree(&out);

  return EXIT_SUCCESS;

}

int test2(char fname[]){
  radiationOutputFields_t out;
  unsigned int i = 0; 
  unsigned int j = 0;
  unsigned int k = 47;
  /* unsigned int kn1; */
  unsigned int t = 0;

  unsigned int offset;

  getAndSetAtmosFieldsFromFile(fname, &out);

  printf("\nTest Results:\n\n");
  /* for(t=0; t<out.ntime; ++t){ */
    /* for(i=0; i<out.nlat; ++i){ */
      /* for(j=0; j<out.nlon; ++j){ */
        /* for(k=0; k<out.npfull; ++k){ */
       /* for(k=out.npfull; k>0; --k){ */
       /*   kn1 = k-1; */
         offset = t*(out.nlat*out.nlon*out.npfull) + i*(out.nlon*out.npfull) + j*(out.npfull) +k;
         printf("DELTAZ(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.DELTAZ[offset]);
         printf("P(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.P[offset]);
         printf("T(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.T[offset]);
         printf("N(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.N[offset]);
         printf("PS(t=%d, lat=%d, lon=%d, pfull=%d )\t=\t%g\n" ,  t, i, j, k, out.PS[offset]);
         
        /* } */
  /*     } */
  /*   } */
  /* } */

  radiationOutputFieldsFree(&out);

  return EXIT_SUCCESS;

}


#ifndef SKIPMAIN
int main(int argc, char* argv[]){
  if(argc != 2){
    fprintf(stderr,"./%s fname.nc\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  if ( /* (test(argv[1])) || */ (test2(argv[1])) ){  
    fprintf(stderr,"Tests FAILED!\n");
    exit(EXIT_FAILURE);
  }
    
  return EXIT_SUCCESS;

}
#endif
      

      
