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

#ifdef MPI_ENABLED
#include <mpi.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <argp.h>

#include "grtcode.h"
#include "TIPS_2011.h"
#include "parseHITRANfile.h"
#include "parseNetcdfRadiation.h"
#include "outputNetcdfSpec.h"
#include "voigt.h"
#include "continuum.h"

/* begin data structure helpers */

/* these happen to be HITRAN_MOLID-1 ...*/
typedef enum MoleculeNumber_t
{
  H2O = 0,
  CO2 = 1,
  O3  = 2,
  N2O = 3,
  CO  = 4,
  CH4 = 5,
  O2  = 6,
  NUM_MOL = 7
} MoleculeNumber_t;

/* end data structures helpers */

/* for argp doc */
const char *argp_program_version = "lbl-dev 0.1";
const char *argp_program_bug_address = "<garrett.wright@noaa.gov>";
static char doc[] = "GFDL style documentation goes >/\n\n\\"
    "<^here.\n\v"
    "Other Documentation goes here.";
/* A description of the arguments we accept. */
#define minNhitfiles 1
#define maxNhitfiles NUM_MOL
static const unsigned int minNargs=minNhitfiles;
static const unsigned int maxNargs=maxNhitfiles;
static char args_doc[] = "-aINPUT.nc -oOUT.nc [molecule concentration specifications] HITFILES";
/* options we'd like to understand */
static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Opens the elevator door"},
  {"quiet",   'q', 0, 0, "Closes the elevator door"},
  {"silent",  's', 0, OPTION_ALIAS },
  {"device",  'd', "VAL", OPTION_ARG_OPTIONAL, "Use gpu implimentation on specifed DEVICE. "
   "\n\tDefault DEVICE is simply GPU0"
   "\n\tIncompatible with --host."
   "\n\t --mpi modifies this flag to prescribe numDevices per node."},
  {"host",    'h', 0, 0, "Use HOST cpu implimentation. \n\t(incompatible with --device)"},
  {"mpi",     'M', 0, 0, "Use MPI: Ranks taken from MPI_Comm_World. (Must be compiled for MPI!)"},
  {"output",  'o', "FILE", 0, "Output to FILE" },
  {"atmos",   'a', "INPUT.NC", 0,"NC file containing model atmosphere."},
  {"minw",    'w', "VAL", OPTION_ARG_OPTIONAL, "minimum Wavenumber (lower bound, inclusive), defaults 1", -3},
  {"maxw",    'W', "VAL", OPTION_ARG_OPTIONAL, "maximum Wavenumber (upper bound, inclusive), defaults 50000", -3},
  {"mint",    't', "VAL", OPTION_ARG_OPTIONAL, "minimum time (lower bound, inclusive), defaults 0", -3},
  {"maxt",    'T', "VAL", OPTION_ARG_OPTIONAL, "maximum time (upper bound, inclusive), defaults 0", -3},
  {"res",     'r', "VAL", OPTION_ARG_OPTIONAL, "Resolution (wavenumber), defaults 1.0", -3},
  {"wings",   'c', "VAL", OPTION_ARG_OPTIONAL, "Wings cutoff (+/- integer wavenumber), defaults 25", -3},
  {"h2o",     '1', "VAL", OPTION_ARG_OPTIONAL, "Water Concentration.  Default reads from INPUT.NC, else supply global PartialPressure(atm) value",-2},
  {"co2",     '2', "VAL", OPTION_ARG_OPTIONAL, "Carbon Dioxide Concentration.  Supply global PartialPressure(atm) value",-2},
  {"o3",      '3', "VAL", OPTION_ARG_OPTIONAL, "Ozone Concentration.  Default reads from INPUT.NC, else supply global PartialPressure(atm) value",-2},
  {"n2o",     '4', "VAL", OPTION_ARG_OPTIONAL, "Nitrous Oxide. Supply global PartialPressure(atm) value",-2},
  {"co",      '5', "VAL", OPTION_ARG_OPTIONAL, "Carbon Monoxide Concentration.  Supply global PartialPressure(atm) value",-2},
  {"ch4",     '6', "VAL", OPTION_ARG_OPTIONAL, "Methane Conentration.  Supply PartialPressure(atm) global value",-2},
  {"o2",      '7', "VAL", OPTION_ARG_OPTIONAL, "Oxygen Concentration.  Supply PartialPressure(atm) global value",-2},
  {"ctm",     'C', 0, 0, "Enables the continuum codes for testing",-1},
  {0}
};

/* our args structure */
struct arguments {
  char *atmos;                  /* input atmos fname */
  char *hitfiles[maxNhitfiles];  /* hit1 hit2...*/
  int nhitfiles;
  int nmolConc;
  int nmolConcOver;
  int silent;
  int verbose;
  int host;
  int device;
  int mpi;
  int t,T,w,W;
  int ctm;
  double res;
  int wingBreadth;
  double h2o,co2,o3,n2o,co,ch4,o2;  
  char *output_file;
};

/* parser */
static double parse_MolecConc(char *arg)
{
  double res;
  if (arg[0] == 'a' )
  {  /* hook to explicitly yield from "nc" file */
    res = -1;
  }
  else if (isalpha(arg[0]))
  {
    fprintf(stderr,"The supplied character (%c) "
            "for overriding a molecule concentration is not understood. Review args, aborting.\n",
            arg[0]);
    exit(EXIT_FAILURE);
  }            
  else
  {
    /* this should probably be changed to a strtod call with err checks later*/
    res = atof(arg);
  }
  return res;
}

static error_t parse_opt( int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our args structure. */
  struct arguments* arguments = (struct arguments*)state->input;
  
  switch(key)
  {
    case 'q':case 's':
      arguments->silent = 1;
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'o':
      arguments->output_file = arg;
      break;
    case 'h':
      arguments->host=1;
      break;
    case 'd':
      arguments->device=atoi(arg);
      break;
    case 'M':
      arguments->mpi=1;
      break;
    case 'a':
      arguments->atmos = arg;
      break;
      /* user prescribed bounds */
    case 'w':
      arguments->w = atoi(arg);
      break;
      /* user prescribed bounds */
    case 'W':
      arguments->W = atoi(arg);
      break;
    case 't':
      arguments->t = atoi(arg);
      break;
      /* user prescribed bounds */
    case 'T':
      arguments->T = atoi(arg);
      break;
    case 'r':
      arguments->res = atof(arg);
      break;      
    case 'c':
      arguments->wingBreadth = atoi(arg);
      break;
      /* Molecular Concentrations */
    case '1':      
      arguments->h2o = parse_MolecConc(arg);
      arguments->nmolConc++;
      break;
    case '2':
      arguments->co2 = parse_MolecConc(arg);
      arguments->nmolConc++;
      if (arguments->co2 == -1)
      {
        fprintf(stderr,"NC read is not currently supported for this molecule (%c), aborting..\n",key);
        exit(EXIT_FAILURE);
      }
      break;
    case '3':
      arguments->o3 = parse_MolecConc(arg);
      arguments->nmolConc++;
      break;
    case '4':
      arguments->n2o = parse_MolecConc(arg);
      arguments->nmolConc++;
      if (arguments->n2o == -1)
      {
        fprintf(stderr,"NC read is not currently supported for this molecule (%c), aborting..\n",key);
        exit(EXIT_FAILURE);
      }
      break;
    case '5':
      arguments->co = parse_MolecConc(arg);
      arguments->nmolConc++;
      if (arguments->co == -1)
      {
        fprintf(stderr,"NC read is not currently supported for this molecule (%c), aborting..\n",key);
        exit(EXIT_FAILURE);
	}
      break;
    case '6':
      arguments->ch4 = parse_MolecConc(arg);
      arguments->nmolConc++;
      if (arguments->ch4 ==-1)
      {
        fprintf(stderr,"NC read is not currently supported for this molecule (%c), aborting..\n",key);
        exit(EXIT_FAILURE);
      }
      break;
    case '7':
      arguments->o2 = parse_MolecConc(arg);
      arguments->nmolConc++;
      if (arguments->o2 == -1)
      {
        fprintf(stderr,"NC read is not currently supported for this molecule (%c), aborting..\n",key);
        exit(EXIT_FAILURE);
      }
      break;
    case 'C':
      arguments->ctm = 1;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= maxNargs )
      {
        fprintf(stderr,"too many arguments can make a parsern go crazy\n");
	  argp_usage(state);
	}
      arguments->hitfiles[state->arg_num] = arg;
      arguments->nhitfiles++;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < minNargs )
	{
	  fprintf(stderr,"too few arguments is just not enough\n");
	  argp_usage( state );
	}
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* /the/ argp pargser */
static struct argp argp = {options, parse_opt, args_doc, doc };

static void setGlobalPartialPres(double val, 
                                 REAL_t* PS,
                                 unsigned int molId,
                                 size_t ntime,
                                 size_t nlat,
                                 size_t nlon,
                                 size_t nlvl)
{
  /* DEPENDS: PS is stored C order PS[time][lats][lons][mols][z]  */
  unsigned int itr;
  unsigned int lat;
  unsigned int lon;
  unsigned int time;
  size_t ps_off;
  const size_t ps_off_mol = (molId-1)*nlvl;

  for(time=0; time<ntime; ++time)
  {
    for(lat=0; lat<nlat; ++lat)
    {
      for(lon=0; lon<nlon; ++lon)
      {
        /* get an offset for this t,lat,lon,mol */
        ps_off = time*nlat*nlon*NUM_MOL*nlvl + lat*nlon*NUM_MOL*nlvl + lon*NUM_MOL*nlvl + ps_off_mol;
        /* itr over levels */
        for(itr=0; itr<nlvl; ++itr)
        {
          PS[ps_off + itr] = val;
        }
      }
    }
  }
}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t idealGasNumberDensity( const REAL_t P_atm,
                              const REAL_t V_cm3,
                              const REAL_t T_k){
  const REAL_t R = 82.057338; /* (cm^3*atm*K^-1*mol^-1) */  
  const REAL_t AV = 6.022E23;
  return P_atm * V_cm3 * AV / (R * T_k );
}

static void setGlobalNumberDensity(REAL_t* const N,    /* molecules/cm^3 */
                                   REAL_t const * const PartialPres,  /* atm */
                                   REAL_t const * const T,  /* kelvin */
                                   const unsigned int hitranMolId,
                                   const size_t ntime,
                                   const size_t nlat,
                                   const size_t nlon,
                                   const size_t nlvl)
{
  /* DEPENDS: N is stored C order N[time][lats][lons][mols][z]  */
  unsigned int itr;
  unsigned int lat;
  unsigned int lon;
  unsigned int time;
  size_t off;
  size_t n_off;
  const size_t n_off_mol = (hitranMolId-1)*nlvl;
  
  for(time=0; time<ntime; ++time)
  {
    for(lat=0; lat<nlat; ++lat)
    {
      for(lon=0; lon<nlon; ++lon)
      {
        /* get an offset for this t,lat,lon,mol */
        off = time*nlat*nlon*nlvl + lat*nlon*nlvl + lon*nlvl;
        n_off = time*nlat*nlon*NUM_MOL*nlvl + lat*nlon*NUM_MOL*nlvl + lon*NUM_MOL*nlvl + n_off_mol;
        /* itr over levels */
        for(itr=0; itr<nlvl; ++itr)
        {
          /* uses ideal gas law, but there is an implied "unit area" so V became 1 */
          /* DEPENDS units */
          N[n_off + itr] = idealGasNumberDensity( PartialPres[n_off+itr], 1.,  T[off+itr]);
        }
      }
    }
  }
}

static void checkMolConfig(struct arguments* args,
                           const unsigned int molid,
                           radiationOutputFields_t* atmosData,
                           const int time)
{
  REAL_t* PS = atmosData->PS;
  const size_t nlat = atmosData->nlat;
  const size_t nlon = atmosData->nlon;
  const size_t nlvl = atmosData->npfull;
  
  
  int abort=0;
  switch(molid)
  {
    case 1:
      if (args->h2o ==0)
      {
        abort=1;
      }
      else if(args->h2o > 0)
      {
        setGlobalPartialPres(args->h2o, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    case 2:
      if (args->co2 ==0)
      {
        abort=1;
      }
      else if(args->co2 > 0)
      {
        setGlobalPartialPres(args->co2, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    case 3:
      if (args->o3 ==0)
      {
        abort=1;
      }
      else if(args->o3 > 0)
      {
        setGlobalPartialPres(args->o3, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    case 4:
      if (args->n2o ==0)
      {
        abort=1;
      }
      break;
    case 5:
      if (args->co ==0)
      {
        abort=1;
      }
      else if(args->co > 0)
      {
        setGlobalPartialPres(args->co, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    case 6:
      if (args->ch4 ==0)
      {
        abort=1;
      }
      else if(args->ch4 > 0)
      {
        setGlobalPartialPres(args->ch4, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    case 7:
      if (args->o2 ==0)
      {
        abort=1;
      }
      else if(args->o2 > 0)
      {
        setGlobalPartialPres(args->o2, PS, molid, time, nlat, nlon, nlvl);
      }
      break;
    default:
      abort=-1;
      break;
  }

  if (abort==-1)
  {
    fprintf(stderr,"This Hitfiles MolId (%d) does not appear to be supported yet.\n", molid);
    exit(EXIT_FAILURE);
  }

  else if (abort!=0)
  {
    fprintf(stderr,"This Hitfile MolId (%d) does not appear to match any of the provided Mol Concentrations."
            " Probably missing a Hitfile, wrong hitfile, or extra Mol Concentration specified.\n", molid);
    exit(EXIT_FAILURE);
  }

  /* else for any molecule that makes it this far we will need it's number density */
  setGlobalNumberDensity(atmosData->N,    /* molecules/cm^3 */
                         PS,  /* atm */
                         atmosData->T,  /* kelvin */
                         molid,
                         time,
                         nlat,
                         nlon,
                         nlvl);

  
  return;
}

#define MAX_NUM_SPECTRAL_LINES 524288  /* 2^19 */
const REAL_t TREF = 296.0;
/* const REAL_t c1 = 1.191042869E-8; /\* W/m2*sr*cm-4 noaa units*\/ /\* used in rad solver *\/ */
const REAL_t c2 = 1.4387686;  /* to match RFM4.3, and Rothman Paper, hc/k in cm K , h:ergs s , c: cm/s , k [ergs/K]*/
/* const REAL_t corK =1.66;      /\* used in rad solver *\/ */
/* for quick reference */
/* const REAL_t c1 = 1.191042869E-16;  /\* first radiation constant W*m2*sr^1*cm^-1 , SI*\/ */
/* const REAL_t c2 = 1.4387770E-2;  /\* second radiation constant meters Kelvin from SI *\/ */

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#ifndef M_1_PI
#define M_1_PI 0.31830988618379067154/* 1/pi */
#endif

#ifndef MAXNSTREAMS
#define MAXNSTREAMS 2
#endif

#ifdef __NVCC__
/* begin gpu helpers */

#undef FORCE_KERNEL_CHECK
/* #define FORCE_KERNEL_CHECK */

#undef EVENTS
/* #define EVENTS */

#include "cudaHelpers.cuh"

#endif
/* end gpu helpers */



/******************** begin absorption calcs *****************************/  

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t getMolarMass(const int hitranMolId){
  REAL_t res;
  switch(hitranMolId)
  {
    case 1:     /* h2o */
      res = 18.01528;  /* g/mol */;
      break;
    case 2:     /* co2 */
      res = 44.01;
      break;
    case 3:  /* o3 */
      res = 48.;
      break;
    case 4:  /* n2o */
      res = 44.013;
      break;
    case 5:  /* co */
      res = 28.01;
      break;
    case 6:  /* ch4 */
      res = 16.04;
      break;
    case 7:  /* o2 */
      res = 32.;
      break;
    default:
/* #if !defined(__CUDA_ARCH__) */
/*       fprintf(stderr,"Error, the molecular with (0-based) molId=%d is not implimented in getMolarMass," */
/*               " something is probably very very wrong. Aborting\n.", hitranMolId); */
/*       exit(EXIT_FAILURE); */
/* #else */
      assert(0);  /* mol not implimented */
/* #endif */
      break;
  }
         
  return res;
}

/* Q found as QT by method of TIPS_2011(.cu) */

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Q(const uint8_t molId,
         const REAL_t T,
         const uint8_t iso)
{
  float gsi;
  REAL_t Qt;
  QT(molId,
     T, /* input temp */
     iso, /* an isotopeIndex */
     &gsi, /* state independent nuclear degeneracyfactor */
     &Qt);       /*  Total Internal Partition Function */
  /* /\* dbg *\/ printf("Q(molId=%d, T=%f, isoIdx=%d) = %f \n", molId,T,iso,Qt); */
  return Qt;
}



#ifdef __NVCC__
__host__ __device__
#endif
REAL_t compute_gamma(const REAL_t P,
                     const REAL_t T,
                     const float Yself,
                     const  float Yair,
                     const  float n,
                     const  REAL_t Ps)
{
  const REAL_t gam = pow( (TREF/T), (REAL_t)n ) * ( ((REAL_t)Yair) *(P-Ps)+((REAL_t)Yself)*Ps);
  /* /\* dbg *\/ printf("compute_gamma(P=%f, T=%f, Yself=%f, Yair=%f, n=%f, Ps=%f) = %f\n", P, T, Yself, Yair, n, Ps, gam); */
  return gam;
}


#ifdef __NVCC__
__host__ __device__
#endif
REAL_t pressureShiftCorrection(const REAL_t Vnn,
                               const float d,
                               const REAL_t P)
{
  return (Vnn + ((REAL_t)d) * P );
}


#ifdef __NVCC__
__host__ __device__
#endif
REAL_t lorentzianKernel( const REAL_t F,
                         const REAL_t gam,
                         const REAL_t gam2,
                         const REAL_t PshiftCorrection)
{
  const REAL_t del = (F - PshiftCorrection);

  return ( ((REAL_t)M_1_PI) * gam / ( gam2 + del*del) );
}


#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Snn_partialCorrection(const uint8_t molId,
                             uint8_t iso,
                             REAL_t Vnn,
                             float En,
                             REAL_t Snn_ref)
{
  return Snn_ref * Q(molId,TREF,iso) /
      ( exp(-c2*En/TREF)  * ( 1 - exp(-c2*(Vnn/TREF) ) ) ) ;
}

#ifdef __NVCC__
__host__ __device__
#endif
REAL_t Snn_Tcorrection(const uint8_t molId,
                       REAL_t T,
                       uint8_t iso,
                       REAL_t Vnn,
                       float En,
                       REAL_t Snn_partial)
{
  return ( Snn_partial / Q(molId,T,iso) ) * ( ((REAL_t)1) - exp(-c2*Vnn/T) ) * exp(-c2*En/T) ;
}


#ifdef __NVCC__
__inline__ __host__ __device__
#endif
REAL_t Knn_Lor(const REAL_t snn,
             const REAL_t F,
             const REAL_t gam,
             const REAL_t gam2,
             const REAL_t PshiftCorrection)
{
  return ( snn * lorentzianKernel(F,gam,gam2,PshiftCorrection) );
}


#ifdef __NVCC__
__inline__ __host__ __device__
#endif
REAL_t tau_Lor(const REAL_t Snn,
           const REAL_t F,
           const REAL_t gam,
           const REAL_t gam2,
           const REAL_t PshiftCorrection,
           const REAL_t u,          /* mol cm-2 */
           const REAL_t pathlength)  /* should match units of u, cm */
{
  return ( (pathlength*u) * Knn_Lor( Snn,F,gam,gam2,PshiftCorrection) ) ;
}

#ifdef __NVCC__
__inline__ __host__ __device__
#endif
REAL_t Knn_Voigt(const REAL_t snn,
                 const REAL_t vnn,
                 const REAL_t F,
                 const REAL_t gam,
                 const REAL_t gam2,
                 const REAL_t PshiftCorrection,
                 const REAL_t eta,
                 const REAL_t alphad)
{
  const float ly = lorentzianKernel(F, gam, gam2, PshiftCorrection);
  const float gy = gauKernel(F, vnn, alphad);
  /* /\* dbg *\/ printf("Lor Kernel: %f \t Gau Kernel(v=%f, vnn=%f, alphad= %f)= %f \n",ly,F,vnn,alphad, gy); */
      
  return ( snn * pseudoVoigt(eta, ly, gy) );
}


#ifdef __NVCC__
__inline__ __host__ __device__
#endif
REAL_t tau_Voigt(const REAL_t Snn,
                  const REAL_t Vnn,
                  const REAL_t F,
                  const REAL_t gam,
                  const REAL_t gam2,
                  const REAL_t PshiftCorrection,
                  const REAL_t eta,
                  const REAL_t alphad,
                  const REAL_t u,          /* mol cm-2 */
                  const REAL_t pathlength)  /* should match units of u, cm */
{  
  return ( (pathlength*u) * Knn_Voigt( Snn,Vnn, F, gam, gam2, PshiftCorrection, eta, alphad) ) ;
}

#ifdef __NVCC__
__global__
void eval_gamma( const unsigned int numLayers,
                 const unsigned int nL,
                 const REAL_t* const P,
                 const REAL_t* const T,
                 const REAL_t* const Ps,
                 const float* const Yself,
                 const float* const Yair,
                 const float* const n,
                 REAL_t* Gam)
{
  unsigned int lyr;
  unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;
  if (ltid < nL)
  {
#pragma unroll
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      Gam[lyr*nL + ltid] = compute_gamma(P[lyr],
                                         T[lyr],
                                         Yself[ltid],
                                         Yair[ltid],
                                         n[ltid],
                                         Ps[lyr]);
    }
  }
  return;
}

__global__
void eval_pShift(const unsigned int numLayers,
                 const unsigned int nL,
                 const REAL_t* const P,
                 const REAL_t* const Vnn,
                 const float* const d,
                 REAL_t* const PShift)
{
  unsigned int lyr;
  unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;
  if (ltid < nL)
  {
#pragma unroll
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      PShift[lyr*nL + ltid] = pressureShiftCorrection(Vnn[ltid],
                                                      d[ltid],
                                                      P[lyr]);
    }
  }
  return;
}

__global__ void pre_eval_Snn(const unsigned int nL,
                             const uint8_t molId,
                             const uint8_t* const iso,
                             const REAL_t* const Vnn,
                             const float* const En,
                             REAL_t* Snn_ref)
{
  int ltid = blockIdx.x*blockDim.x + threadIdx.x;
  if (ltid < nL)
  {
    Snn_ref[ltid] = Snn_partialCorrection(molId,
                                          iso[ltid],
                                          Vnn[ltid],
                                          En[ltid],
                                          Snn_ref[ltid]);
  }
  return;
}

__global__ void eval_Snn_correction(const unsigned int numLayers,
                                    const unsigned int nL,
                                    const uint8_t molId,
                                    const REAL_t* const T,
                                    const uint8_t* const iso,
                                    const REAL_t* const Vnn,
                                    const float* const En,
                                    const REAL_t* const Snn_partial,
                                    REAL_t* const S)
{
  unsigned int lyr;
  unsigned int ltid = blockIdx.x*blockDim.x + threadIdx.x;
  if (ltid < nL)
  {
#pragma unroll
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      S[lyr*nL + ltid] = Snn_Tcorrection(molId,
                                         T[lyr],
                                         iso[ltid],
                                         Vnn[ltid],
                                         En[ltid],
                                         Snn_partial[ltid]);
    }
  }
  return;
}


__global__ void eval_profile(const unsigned int molId,  /* hitran, 1 based */
                             const unsigned int nL,
                             const int nF,
                             const REAL_t loWn,
                             const REAL_t resolution,
                             const unsigned int numLayers,
                             const unsigned int breadth,

                             const REAL_t* const T,
                             REAL_t* Vnn,
                             
                             REAL_t* Gam,
                             REAL_t* PShift,
                             REAL_t* S,
                             REAL_t* tauU_d,
                             REAL_t* pathlength_d,
                             REAL_t* out)
{
  unsigned int ltid = blockIdx.x * blockDim.x + threadIdx.x;
  if ( ltid<nL )
  {
    int ftid;
    int fcenterid;
    REAL_t f;
    unsigned int lyr;
    unsigned int loffset;
    
    REAL_t gam;
    REAL_t gam2;
    REAL_t pShift;
    REAL_t snn;
    REAL_t tauu;
    REAL_t len;

    const REAL_t molarMass = getMolarMass(molId);  /* g/mol */
    REAL_t temp;
    REAL_t etav;
    REAL_t alphad;
    REAL_t gaufwhm;


    const int fsteps = ceil((REAL_t)breadth/resolution);
    
    /* Find index of nearest frequency to line */
    const REAL_t thisLine = Vnn[ltid] ;
    fcenterid = ( (2*(thisLine-loWn)/resolution) + 1 )/2;
    if (fcenterid<nF){
      
#pragma unroll
      for( lyr=0; lyr< numLayers ; ++lyr ) /*  spatial points */
      {
        loffset = lyr*nL + ltid;
        /* lookup Gamma, pressureShift, Snn at this spatial point for this line*/
        gam = Gam[loffset];
        gam2 = gam*gam;
        pShift = PShift[loffset];
        snn = S[loffset];
        tauu = tauU_d[lyr];  /* zero base the molId */
        len = pathlength_d[lyr];

        temp = T[lyr];
        gaufwhm = gauFWHM(temp,molarMass,thisLine);
        /* /\* dbg *\/ printf("lorfwhm = %f \t gaufwhm = %f \n",gam, gaufwhm); */
        etav = eta(gam, gaufwhm );
        alphad = gauAlphad(temp,molarMass,thisLine);

        
#pragma unroll
        for( ftid= fcenterid - ((int)fsteps); ftid <= fcenterid ; ++ftid)
        {
          if (ftid>=0){
            f = ((REAL_t)ftid)*resolution + loWn;
            /* must use atomics for now, incorrect (race on load-alter-write out[ftid] ) without presorting lines (which we can do later) :) */
            atomicAdd( &(out[lyr*nF + ftid]) , tau_Voigt(snn, thisLine, f, gam, gam2, pShift, etav, alphad, tauu, len) );
          }
        }
#pragma unroll
        for( ftid = fcenterid + ((int)fsteps) ; ftid > fcenterid; --ftid)
        {
          if( ftid<nF){
            f = ((REAL_t)ftid)*resolution + loWn;
            /* must use atomics for now, incorrect (race on load-alter-write out[ftid] ) without presorting lines (which we can do later) :) */
            atomicAdd( &(out[lyr*nF + ftid]) , tau_Voigt(snn, thisLine, f, gam, gam2, pShift, etav, alphad, tauu, len) );
          }
        }

      }
    }
  }
  return;
}
#endif

/* HOST */

void eval_gamma_h(const unsigned int numLayers,
                  const unsigned int nL,
                  const REAL_t* const P,
                  const REAL_t* const T,
                  const REAL_t* const Ps,
                  const float* const Yself,
                  const float* const Yair,
                  const float* const n,
                  REAL_t* Gam)
{
  unsigned int lyr;
  unsigned int ltid;
  for(ltid=0; ltid < nL; ++ltid)
  {
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      /* /\* dbg *\/ printf("eval_gamma_h( lyr= %d, ltid=%d) \n",lyr,ltid); */
      Gam[lyr*nL + ltid] = compute_gamma(P[lyr],
                                         T[lyr],
                                         Yself[ltid],
                                         Yair[ltid],
                                         n[ltid],
                                         Ps[lyr]);
    }
  }
  return;
}

void eval_pShift_h(const unsigned int numLayers,
                  const unsigned int nL,
                  const REAL_t* const P,
                  const REAL_t* const Vnn,
                  const float* const d,
                  REAL_t* const PShift)
{
  unsigned int lyr;
  unsigned int ltid;
  for(ltid=0; ltid < nL ; ++ltid)
  {
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      PShift[lyr*nL + ltid] = pressureShiftCorrection(Vnn[ltid],
                                                      d[ltid],
                                                      P[lyr]);
    }
  }
  return;
}

void pre_eval_Snn_h(const unsigned int nL,
                    const uint8_t molId,
                    const uint8_t* const iso,
                    const REAL_t* const Vnn,
                    const float* const En,
                    REAL_t* Snn_ref)
{
  unsigned int ltid;
  for(ltid=0; ltid < nL ; ++ltid)
  {
    Snn_ref[ltid] = Snn_partialCorrection(molId,
                                          iso[ltid],
                                          Vnn[ltid],
                                          En[ltid],
                                          Snn_ref[ltid]);
  }
  return;
}

void eval_Snn_correction_h(const unsigned int numLayers,
                           const unsigned int nL,
                           const uint8_t molId,
                           const REAL_t* const T,
                           const uint8_t* const iso,
                           const REAL_t* const Vnn,
                           const float* const En,
                           const REAL_t* const Snn_partial,
                           REAL_t* const S)
{
  unsigned int lyr;
  unsigned int ltid;
  for(ltid=0; ltid < nL ; ++ltid )
  {
    for( lyr = 0 ; lyr < numLayers; ++lyr )
    {
      S[lyr*nL + ltid] = Snn_Tcorrection(molId,
                                         T[lyr],
                                         iso[ltid],
                                         Vnn[ltid],
                                         En[ltid],
                                         Snn_partial[ltid]);
    }
  }
  return;
}


void eval_profile_h(const unsigned int molId,
                    const unsigned int nL,
                    const int nF,
                    const REAL_t loWn,
                    const REAL_t resolution,
                    const unsigned int numLayers,
                    const unsigned int breadth,

                    const REAL_t* const T,
                    REAL_t* Vnn,

                    REAL_t* Gam,
                    REAL_t* PShift,
                    REAL_t* S,
                    REAL_t* tauU,
                    REAL_t* pathlength,
                    REAL_t* out)
{
  unsigned int ltid;
  for(ltid=0; ltid<nL ; ++ltid)
  {
    int ftid;
    int fcenterid;
    REAL_t f;
    unsigned int lyr;
    unsigned int loffset;
    
    REAL_t gam;
    REAL_t gam2;
    REAL_t pShift;
    REAL_t snn;
    REAL_t tauu;
    REAL_t len;

    const REAL_t molarMass = getMolarMass(molId);
    REAL_t temp;
    REAL_t etav;
    REAL_t alphad;
    REAL_t gaufwhm;
    
    const int fsteps = ceil((REAL_t)breadth/resolution);
    
    /* Find index of nearest frequency to line */
    const REAL_t thisLine = Vnn[ltid] ;
    fcenterid = ( (2*(thisLine-loWn)/resolution) + 1 )/2;
    if (fcenterid<nF)
    {      
      for( lyr=0; lyr< numLayers ; ++lyr ) /*  spatial points */
      {
        loffset = lyr*nL + ltid;
        /* lookup Gamma, pressureShift, Snn at this spatial point for this line*/
        gam = Gam[loffset];
        gam2 = gam*gam;
        pShift = PShift[loffset];
        snn = S[loffset];
        tauu = tauU[lyr];
        len = pathlength[lyr];

        temp = T[lyr];
        gaufwhm = gauFWHM(temp,molarMass,thisLine);
        
        /* /\* dbg *\/ printf("lorfwhm = %f \t gaufwhm(temp=%f, molarMass=%f, v0=%f) = %f \n",gam, temp, molarMass, thisLine, gaufwhm); */
        etav = eta(gam, gaufwhm );
        alphad = gauAlphad(temp,molarMass,thisLine);
        
        for( ftid= fcenterid - ((int)fsteps); ftid <= fcenterid ; ++ftid)
        {
          if (ftid>=0){
            f = ((REAL_t)ftid)*resolution + loWn;
            out[lyr*nF + ftid] += tau_Voigt(snn, thisLine, f, gam, gam2, pShift, etav, alphad, tauu, len);
          }
        }

        for( ftid = fcenterid + ((int)fsteps) ; ftid > fcenterid; --ftid)
        {
          if( ftid<nF){
            f = ((REAL_t)ftid)*resolution + loWn;
            out[lyr*nF + ftid] += tau_Voigt(snn, thisLine, f, gam, gam2, pShift, etav, alphad, tauu, len);
          }
        }
        
        /* /\* dbg *\/ printf("eval_profile(lyr=%d, ltid= %d, gam=%f, S=%e, tauu=%f ): centerF=%f centerV = %f \n", lyr, ltid, gam, snn, tauu, f, out[lyr*nF + fcenterid]); */
        
      }
    }
  }
  return;
}

int host_optics_perMol(const uint8_t molId,
                       const unsigned int nL,
                       const REAL_t loWn,
                       /* const REAL_t hiWn, */
                       const unsigned int nF,
                       const REAL_t resolution,
                       const int breadth,
                       const unsigned int numLayers,
                       REAL_t* out_h,
                       const REAL_t* const T_h,
                       const REAL_t* const P_h,
                       const uint8_t* const iso,
                       const REAL_t* const Vnn,
                       const REAL_t* const Snn_ref,
                       const float* const Yair,
                       const float* const Yself,
                       const float* const En,
                       const float* const n,
                       const float* const d,
                       REAL_t* const TauU_h,   /* number density precomputed in atmosData */
                       REAL_t* const pathLength_h,  /* pathlength precomputed in atmosData */
                       const REAL_t* const PS_h,   /* partial pres precomputed in atmosData */
                       RefLinePtrs_t* const Lines_h,
                       OpticsBufPtrs_t* const OptBuf_h)
{
  
  /* Initilization */
  initTIPS();

  /* declares */
  printf( "HITRAN molId:=%d\n", molId);
  
  /* const unsigned int nL=Lines_h->nLines; */
      
  uint8_t* iso_h=Lines_h->iso;
  REAL_t* Vnn_h=Lines_h->Vnn;
  REAL_t* Snn_ref_h=Lines_h->Snn_ref;
  float* Yair_h=Lines_h->Yair;
  float* Yself_h=Lines_h->Yself;
  float* En_h=Lines_h->En;
  float* n_h=Lines_h->n;
  float* d_h=Lines_h->d;

  /* per mol */
  memcpy(iso_h, iso, nL*sizeof(uint8_t) );
  memcpy(Vnn_h, Vnn, nL*sizeof(REAL_t) );
  memcpy(Snn_ref_h, Snn_ref, nL*sizeof(REAL_t) );
  memcpy(Yair_h, Yair, nL*sizeof(float) );
  memcpy(Yself_h, Yself, nL*sizeof(float) );
  memcpy(En_h, En, nL*sizeof(float) );
  memcpy(n_h, n, nL*sizeof(float) );
  memcpy(d_h, d, nL*sizeof(float) );

  
  /* precompute/staging */
  REAL_t* Gam_h=OptBuf_h->Gam;
  REAL_t* PShift_h=OptBuf_h->PShift;
  REAL_t* S_h=OptBuf_h->S;

  /*  execute the pre_eval_Snn kernel */
  printf("%s Kernel Launch.. %d lines ... ","pre_eval_Snn_h", nL);
  pre_eval_Snn_h( nL, molId, iso_h, Vnn_h, En_h, Snn_ref_h); /* replaces Snn_ref_h! */
  printf("..done!\n");

  /*  execute the gamma kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_gamma_h", nL, numLayers);
  eval_gamma_h( numLayers, nL, P_h, T_h, PS_h, Yself_h, Yair_h, n_h, Gam_h);
  printf("..done!\n");

  /*  execute the pshift kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_pShift_h", nL, numLayers);
  eval_pShift_h( numLayers, nL, P_h, Vnn_h, d_h, PShift_h );
  printf("..done!\n");

  /*  execute the Snn kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_Snn_correction_h", nL,numLayers);
  eval_Snn_correction_h( numLayers, nL, molId, T_h, iso_h, Vnn_h, En_h, Snn_ref_h, S_h);
  printf("..done!\n");

  /* /\* dbg *\/ eval_Snn_h(numLayers, nL, molId, T_h, iso_h, Vnn_h, En_h, Snn_ref_h, S_h); */
  
  /*  execute the profile kernel */
  printf("%s Kernel Launch.. %d lines, %d fsamples ... ","eval_profile_h", nL, nF);
  eval_profile_h( molId, nL, nF, loWn, resolution, numLayers, breadth, T_h, Vnn_h, Gam_h, PShift_h, S_h, pathLength_h, TauU_h, out_h);
  printf("..done!\n");
  
  /* end per Mol computations*/
  return EXIT_SUCCESS;

}

int host_optics_init(unsigned int numLayers,
                     unsigned int numBufLines,
                     RefLinePtrs_t* L_h,
                     OpticsBufPtrs_t* OpticBuf_h,
                     RefLine_flags_t flags)
{
  /* lines scratch */
  *L_h = allocHost(numBufLines,flags);

  /* Optics Buf  */
  OpticBuf_h->Gam = (REAL_t*)malloc( numLayers*numBufLines*sizeof(REAL_t) );
  OpticBuf_h->PShift = (REAL_t*)malloc( numLayers*numBufLines*sizeof(REAL_t) );
  OpticBuf_h->S = (REAL_t*)malloc( numLayers*numBufLines*sizeof(REAL_t) );

  return EXIT_SUCCESS;
}

int host_optics_free(RefLinePtrs_t* L_h,
                     OpticsBufPtrs_t* OpticBuf_h,
                     RefLine_flags_t flags)
{

  /* Lines */
  freeHost(*L_h, flags);
  
  /* precompute/staging */
  free(  OpticBuf_h->Gam ) ;
  free(  OpticBuf_h->PShift ) ;
  free(  OpticBuf_h->S ) ;

  return EXIT_SUCCESS;

}


int host_launch(const unsigned int numMols,
                /* char**  molFnames, */
                RefLinePtrs_t L[],
                const REAL_t loWn,
                /* const REAL_t hiWn, */
                const unsigned int nF,
                const REAL_t resolution,
                const unsigned int wingBreadth,
                radiationOutputFields_t* atmosData,
                REAL_t* const out,
                int time,
                int lat,
                int lon){

  unsigned int mol;
  unsigned int numLayers=atmosData->npfull;

  /* sanity check */
  assert(numMols<=NUM_MOL);

  /* shift all the atmos data over by time, lat,lon,layers */
  const size_t idx = (time*atmosData->nlon*atmosData->nlat + lat*atmosData->nlon + lon) * atmosData->npfull ;
  REAL_t* T = &(atmosData->T[idx]);
  REAL_t* P = &(atmosData->P[idx]);
  REAL_t* DELTAZ = &(atmosData->DELTAZ[idx]);
  REAL_t* N = &(atmosData->N[ idx*NUM_MOL ]);
  REAL_t* PS = &(atmosData->PS[ idx*NUM_MOL ]);

  OpticsBufPtrs_t OpticsBuf_h;
  RefLinePtrs_t LinesBuf_h;
  
  RefLine_flags_t flags= {((unsigned int) -1),1,0}; /* host cuda malloc default, host=True, device=false */
  
  host_optics_init(numLayers, MAX_NUM_SPECTRAL_LINES, &LinesBuf_h, &OpticsBuf_h,flags);

  for(mol=0; mol < numMols; ++mol ){
    
    host_optics_perMol(L[mol].mol,
                       L[mol].nLines,
                       loWn,
                       nF,
                       resolution,
                       wingBreadth,
                       numLayers,
                       out,
                       T,
                       P,
                       L[mol].iso,
                       L[mol].Vnn,
                       L[mol].Snn_ref,
                       L[mol].Yair,
                       L[mol].Yself,
                       L[mol].En,
                       L[mol].n,
                       L[mol].d,
                       &(N[ (L[mol].mol-1) * atmosData->npfull ]),
                       DELTAZ,
                       /* PS, */
                       &(PS[ (L[mol].mol-1) * atmosData->npfull ]), /* offest [mol] into idx
                                                                     * so PS should be t,lat,lon,mol,z*/
                       &(LinesBuf_h),
                       &(OpticsBuf_h) ); 
  }
  
  /* cleanup all, this is dirty,  */
  host_optics_free(&LinesBuf_h, &OpticsBuf_h, flags);

  return EXIT_SUCCESS;
}

/* END HOST */

#ifdef __NVCC__

__host__ int device_atmos_init(const unsigned int numLayers,
                               const unsigned int numMols,
                               const REAL_t* const T_h,
                               const REAL_t* const P_h,
                               const REAL_t* const N_h,
                               const REAL_t* const Z_h,
                               const REAL_t* const Ps_h,                               
                               REAL_t** T_d,
                               REAL_t** P_d,
                               REAL_t** N_d,
                               REAL_t** Z_d,
                               REAL_t** Ps_d)
{

  printf("\nInitializing atmosphere on device:\n");

  /* mallocs */
  printf("\tMallocing atmospheric data on device..");
  HANDLE_ERROR( cudaMalloc(T_d, (numLayers)*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(P_d, numLayers*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(N_d, numMols*numLayers*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(Z_d, numLayers*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(Ps_d, (numLayers+1)*numMols*sizeof(REAL_t) ) );
  printf(".done!\n");

  /* memcpy 2 dev */
  printf("\tCopying atmospheric data to device..\n");
  printf("\tMemcpy..T_h -> T_d");
  HANDLE_ERROR( cudaMemcpy(*T_d, T_h, numLayers*sizeof(REAL_t), cudaMemcpyHostToDevice ) ); /* +=1 stores surface T at last place in array */
  printf(".done!\n");
  printf("\tMemcpy..P_h -> P_d");
  HANDLE_ERROR( cudaMemcpy(*P_d, P_h, numLayers*sizeof(REAL_t), cudaMemcpyHostToDevice ) );
  printf(".done!\n");
  printf("\tMemcpy.. N_h -> N_d");
  HANDLE_ERROR( cudaMemcpy(*N_d, N_h, numMols*numLayers*sizeof(REAL_t), cudaMemcpyHostToDevice ) );
  printf(".done!\n");
  printf("\tMemcpy.. Z_h -> Z_d");
  HANDLE_ERROR( cudaMemcpy(*Z_d, Z_h, numLayers*sizeof(REAL_t), cudaMemcpyHostToDevice ) );
  printf(".done!\n");
  printf("\tMemcpy.. Ps_h -> Ps_d");
  HANDLE_ERROR( cudaMemcpy(*Ps_d, Ps_h, numLayers*NUM_MOL*sizeof(REAL_t), cudaMemcpyHostToDevice ) );
  printf(".done!\n");

  printf("Initializing atmosphere on device successful.\n\n");
 
  return EXIT_SUCCESS;
}

__host__ int device_atmos_free(REAL_t* T_d,
                               REAL_t* P_d,
                               REAL_t* N_d,
                               REAL_t* Ps_d)
{
 
  HANDLE_ERROR(cudaFree(T_d));
  HANDLE_ERROR(cudaFree(P_d));
  HANDLE_ERROR(cudaFree(N_d));
  HANDLE_ERROR(cudaFree(Ps_d));

  return EXIT_SUCCESS;
}


__host__ int device_optics_init(unsigned int numLayers,
                                unsigned int numBufLines,
                                RefLinePtrs_t* L_d,
                                OpticsBufPtrs_t* OpticBuf_d)
{
#ifdef EVENTS
  cudaEvent_t start;
  cudaEvent_t stop;
  float elapsed;
  HANDLE_ERROR( cudaEventCreate(&start) );
  HANDLE_ERROR( cudaEventCreate(&stop) );
  HANDLE_ERROR( cudaEventRecord(start) ); 
#endif

  /* dev lines */
  *L_d = allocDevice(numBufLines);
  
  /* Optics Buf  */
  HANDLE_ERROR( cudaMalloc(&(OpticBuf_d->Gam), numLayers*numBufLines*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(&(OpticBuf_d->PShift), numLayers*numBufLines*sizeof(REAL_t) ) );
  HANDLE_ERROR( cudaMalloc(&(OpticBuf_d->S), numLayers*numBufLines*sizeof(REAL_t) ) );

#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time to device_optics_init: %3.1f ms\n", elapsed);
#endif

#ifdef EVENTS
  /* cleanup event objects */
  HANDLE_ERROR( cudaEventDestroy(start) );
  HANDLE_ERROR( cudaEventDestroy(stop) );
#endif
  
  return EXIT_SUCCESS;
}

__host__ int device_optics_free(RefLinePtrs_t* L_d,
                                OpticsBufPtrs_t* OpticBuf_d)
{

#ifdef EVENTS
  cudaEvent_t start;
  cudaEvent_t stop;
  float elapsed;
  HANDLE_ERROR( cudaEventCreate(&start) );
  HANDLE_ERROR( cudaEventCreate(&stop) );
  HANDLE_ERROR( cudaEventRecord(start) ); 
#endif
  
  /* cleanup */

  /* Lines */
  HANDLE_ERROR( cudaFree( L_d->iso ));
  HANDLE_ERROR( cudaFree( L_d->Vnn ));
  HANDLE_ERROR( cudaFree( L_d->Snn_ref ));
  HANDLE_ERROR( cudaFree( L_d->Yair ));
  HANDLE_ERROR( cudaFree( L_d->Yself ));
  HANDLE_ERROR( cudaFree( L_d->En ));
  HANDLE_ERROR( cudaFree( L_d->n ));
  HANDLE_ERROR( cudaFree( L_d->d ));

  /* precompute/staging */
  HANDLE_ERROR( cudaFree( OpticBuf_d->Gam ));
  HANDLE_ERROR( cudaFree( OpticBuf_d->PShift ));
  HANDLE_ERROR( cudaFree( OpticBuf_d->S ));
  
#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time to device_optics_free: %3.1f ms\n", elapsed);
#endif

#ifdef EVENTS
  /* cleanup event objects */
  HANDLE_ERROR( cudaEventDestroy(start) );
  HANDLE_ERROR( cudaEventDestroy(stop) );
#endif

  return EXIT_SUCCESS;

}


__host__ int device_optics_perMol(cudaStream_t stream,
                                  const uint8_t molId,
                                  const unsigned int nL,
                                  const REAL_t loWn,
                                  /* const REAL_t hiWn, */
                                  const unsigned int nF,
                                  const REAL_t resolution,
                                  const int breadth,
                                  const unsigned int numLayers,
                                  REAL_t* out_d,
                                  const REAL_t* const T_d,
                                  const REAL_t* const P_d,
                                  const uint8_t* const iso,
                                  const REAL_t* const Vnn,
                                  const REAL_t* const Snn_ref,
                                  const float* const Yair,
                                  const float* const Yself,
                                  const float* const En,
                                  const float* const n,
                                  const float* const d,
                                  REAL_t* const TauU_d,
                                  REAL_t* const pathLength_d,
                                  const REAL_t* const PS_d,
                                  RefLinePtrs_t* const Lines_d,
                                  OpticsBufPtrs_t* const OptBuf_d)
{
  
  /* device mem 
   * declare */

  uint8_t* iso_d=Lines_d->iso;
  REAL_t* Vnn_d=Lines_d->Vnn;
  REAL_t* Snn_ref_d=Lines_d->Snn_ref;
  float* Yair_d=Lines_d->Yair;
  float* Yself_d=Lines_d->Yself;
  float* En_d=Lines_d->En;
  float* n_d=Lines_d->n;
  float* d_d=Lines_d->d;

  /* precompute/staging */
  REAL_t* Gam_d=OptBuf_d->Gam;
  REAL_t* PShift_d=OptBuf_d->PShift;
  REAL_t* S_d=OptBuf_d->S;

  /*dbg*/ printf( "Hitran molId:=%d\n", molId);

  /* api setup (easier across different gpus), taken from documentation */
  int dimBlock=0;   /*  The launch configurator returned block size */
  int minGridSize=0; /*  The minimum grid size needed to achieve the maximum occupancy for a full device launch */
  int dimGrid=0;    /*  The actual grid size needed, based on input size */

#ifdef EVENTS
  cudaEvent_t start;
  cudaEvent_t stop;
  float elapsed;
  HANDLE_ERROR( cudaEventCreate(&start) );
  HANDLE_ERROR( cudaEventCreate(&stop) );
  HANDLE_ERROR( cudaEventRecord(start,stream) ); 
#endif

  /* per mol */
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(iso_d, iso, nL*sizeof(uint8_t), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(Vnn_d, Vnn, nL*sizeof(REAL_t), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(Snn_ref_d, Snn_ref, nL*sizeof(REAL_t), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(Yair_d, Yair, nL*sizeof(float), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(Yself_d, Yself, nL*sizeof(float), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(En_d, En, nL*sizeof(float), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(n_d, n, nL*sizeof(float), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");
  printf("Memcpy..");
  HANDLE_ERROR( cudaMemcpyAsync(d_d, d, nL*sizeof(float), cudaMemcpyHostToDevice, stream ) );
  printf(".done!\n");

#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time to cudaMemcpy: %3.1f ms\n", elapsed);
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif


  /* cuda API may produce warning that can be safely ignored depending on sdk version and -W flags*/
  HANDLE_ERROR( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &dimBlock, pre_eval_Snn, 0, ((int)nL) ) );
  /*  Round up according to array size */
  dimGrid = (((int)nL) + dimBlock - 1) / dimBlock;

#ifdef EVENTS
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif
  /*  execute the gamma kernel */
  printf("%s Kernel Launch.. %d lines ... ","pre_eval_Snn", nL);
  printf("using dimBlock = %d and dimGrid = %d  ... ",dimBlock,dimGrid);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  pre_eval_Snn<<< ((unsigned int)dimGrid) , ((unsigned int)dimBlock) , 0 , stream >>>( nL, molId, iso_d, Vnn_d, En_d, Snn_ref_d); /* replaces Snn_ref_d! */
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaPeekAtLastError() );
  printf("..peek-OK..");
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  printf("..done!\n");

#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time in Kernel: %3.1f ms\n", elapsed);
#endif

  /* cuda API will produce warning that can be safely ignored */
  HANDLE_ERROR( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &dimBlock, eval_gamma, 0, ((int)nL) ) );
  /*  Round up according to array size */
  dimGrid = (((int)nL) + dimBlock - 1) / dimBlock;
  
#ifdef EVENTS
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif
  /*  execute the gamma kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_gamma", nL, numLayers);
  printf("using dimBlock = %d and dimGrid = %d  ... ",dimBlock,dimGrid);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  eval_gamma<<< ((unsigned int)dimGrid), ((unsigned int)dimBlock) , 0 , stream >>>( numLayers, nL, P_d, T_d, PS_d, Yself_d, Yair_d, n_d, Gam_d);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaPeekAtLastError() );
  printf("..peek-OK..");
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  printf("..done!\n");
  
#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time in Kernel: %3.1f ms\n", elapsed);
#endif

  /* api setup declares (easier across different gpus), taken from documentation */
  HANDLE_ERROR( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &dimBlock, eval_pShift, 0, 0) );
  /*  Round up according to array size */
  dimGrid = (((int)nL) + dimBlock - 1) / dimBlock;

#ifdef EVENTS
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif

  /*  execute the pshift kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_pShift", nL, numLayers);
  printf("using dimBlock = %d and dimGrid = %d  ... ",dimBlock,dimGrid);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  eval_pShift<<< ((unsigned int)dimGrid) , ((unsigned int)dimBlock),0,stream >>>( numLayers, nL, P_d, Vnn_d, d_d, PShift_d );
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaPeekAtLastError() );
  printf("..peek-OK..");
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  printf("..done!\n");

#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time in Kernel: %3.1f ms\n", elapsed);
#endif

  /* api setup declares (easier across different gpus), taken from documentation */
  HANDLE_ERROR( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &dimBlock, eval_Snn_correction, 0, 0) );
  /*  Round up according to array size */
  dimGrid = (((int)nL) + dimBlock - 1) / dimBlock;

#ifdef EVENTS
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif
  /*  execute the Snn kernel */
  printf("%s Kernel Launch.. %d lines, %d spatial points ... ","eval_Snn_correction", nL,numLayers);
  printf("using dimBlock = %d and dimGrid = %d  ... ",dimBlock,dimGrid);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  eval_Snn_correction<<< ((unsigned int)dimGrid), ((unsigned int)dimBlock),0,stream >>>( numLayers, nL, molId, T_d, iso_d, Vnn_d, En_d, Snn_ref_d, S_d);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaPeekAtLastError() );
  printf("..peek-OK..");
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  printf("..done!\n");

#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time in Kernel: %3.1f ms\n", elapsed);
#endif

  /* api setup declares (easier across different gpus), taken from documentation */
  HANDLE_ERROR( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &dimBlock, eval_profile, 0, 0) );
  /*  Round up according to array size */
  dimGrid = ( ((int)nL) + dimBlock - 1) / dimBlock;

#ifdef EVENTS
  /* kernel */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif
  /*  execute the profile kernel */
  printf("%s Kernel Launch.. %d lines, %d fsamples ... ","eval_profile", nL, nF);
  printf("using dimBlock = %d and dimGrid = %d  ... ",dimBlock,dimGrid);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  eval_profile<<< ((unsigned int)dimGrid) , ((unsigned int)dimBlock),0,stream >>>( molId, nL, nF, loWn, resolution, numLayers, breadth, T_d, Vnn_d, Gam_d, PShift_d, S_d, pathLength_d, TauU_d, out_d);
#ifdef FORCE_KERNEL_CHECK
  HANDLE_ERROR( cudaPeekAtLastError() );
  printf("..peek-OK..");
  HANDLE_ERROR( cudaDeviceSynchronize() );
#endif
  printf("..done!\n");
  
#ifdef EVENTS
  HANDLE_ERROR( cudaEventRecord(stop,stream) );
  HANDLE_ERROR( cudaEventSynchronize(stop));
  HANDLE_ERROR( cudaEventElapsedTime(&elapsed,start,stop) );
  printf( "Time in Kernel: %3.1f ms\n", elapsed);
  /*  */
  HANDLE_ERROR( cudaEventRecord(start,stream) );
#endif
  
  /* end per Mol computations*/

#ifdef EVENTS
  /* cleanup event objects */
  HANDLE_ERROR( cudaEventDestroy(start) );
  HANDLE_ERROR( cudaEventDestroy(stop) );
#endif

  return EXIT_SUCCESS;

}

__host__ static void initStreams(const unsigned int numMols, cudaStream_t** streams, int* nStreams){
  /* CUDA streams */
  int s;

  if (*nStreams ==-1)
  {
    if (*streams==NULL){
      cudaDeviceProp p = checkDeviceProps();
      if (p.concurrentKernels){
        *nStreams = ( numMols <= MAXNSTREAMS ? numMols : MAXNSTREAMS );
      }
      else{
        *nStreams=1;
      }
      fprintf(stderr,"\nUsing %d streams with %d molecules.\n", *nStreams, numMols);
      
      (*streams) =(cudaStream_t*)malloc(sizeof(cudaStream_t)*(*nStreams));
      
      if (*nStreams>1)
      {
        /* create non default streams */
        for(s=0 ; s<(*nStreams) ; ++s){
          HANDLE_ERROR( cudaStreamCreate( &( (*streams)[s]) ) );
        }
      }
      else
      {  /* default */
        (*streams)[0]=NULL;
      }
    }
  }
  /* else pass */
}
  


__host__ int device_launch(int* nStreams,
                           cudaStream_t** streams,
                           const unsigned int numMols,
                           RefLinePtrs_t L[],
                           const REAL_t loWn,
                           /* const REAL_t hiWn, */
                           const unsigned int nF,
                           const REAL_t resolution,
                           const unsigned int wingBreadth,
                           radiationOutputFields_t* atmosData,
                           REAL_t* const out,
                           int time,
                           int lat,
                           int lon){

  unsigned int m;
  unsigned int mol;
  unsigned int numLayers=atmosData->npfull;

  /* Initilization */
  initTIPS_d();
  initStreams(numMols, streams, nStreams);
  
  /* per stream vars */
  int s;
  /* char* hitFname[nStreams]; */
  
  /* sanity check */
  assert(numMols<=NUM_MOL);
  
  /* resolve host pointers */
  /* shift all the atmosData to lat,lon,layer */
  const size_t idx = (time*atmosData->nlon*atmosData->nlat + lat*atmosData->nlon + lon) * atmosData->npfull ;
  REAL_t* T = &(atmosData->T[idx]);
  REAL_t* P = &(atmosData->P[idx]);
  REAL_t* DELTAZ = &(atmosData->DELTAZ[idx]);
  REAL_t* N = &(atmosData->N[idx*NUM_MOL]);
  REAL_t* PS = &(atmosData->PS[idx*NUM_MOL]);
  /* declare device pointers */
  REAL_t* T_d;
  REAL_t* P_d;
  REAL_t* N_d;
  REAL_t* Z_d;
  REAL_t* PS_d;
  REAL_t* out_d;

  /* RefLinePtrs_t L[numMols]; */
  /* init atmos on device */  
  device_atmos_init(numLayers,
                    NUM_MOL,  /* should be the hardcoded(7) */
                    T,
                    P,
                    N,
                    DELTAZ,
                    PS,
                    &T_d,
                    &P_d,
                    &N_d,
                    &Z_d,
                    &PS_d);
                    
  /* out_h should be Zeros or you're gonna have a bad time! */
  HANDLE_ERROR( cudaMalloc(&out_d, (numLayers+1)*nF*sizeof(REAL_t) ) ); /* re-used */
  printf("Memcpy.. out_h -> out_d");
  HANDLE_ERROR( cudaMemcpy(out_d, out, numLayers*nF*sizeof(REAL_t), cudaMemcpyHostToDevice ) );
  printf(".done!\n");

  /* setup buffers for re-use under streaming */
  OpticsBufPtrs_t OpticsBuf_d[*nStreams];
  RefLinePtrs_t LinesBuf_d[*nStreams];
  for (s=0; s<*nStreams; ++s){
    device_optics_init(numLayers,
                       MAX_NUM_SPECTRAL_LINES,
                       &(LinesBuf_d[s]),
                       &(OpticsBuf_d[s]) );
  }
  
  for (m=0; m < (numMols+(*nStreams-1)); m+=*nStreams )
  {    
    /* perform LBL optical calcs */
    for ( s=0; s<*nStreams; ++s){
      mol=m+s;
      if ( mol < numMols){
        device_optics_perMol((*streams)[s],
                             L[mol].mol,
                             L[mol].nLines,
                             loWn,
                             /* hiWn, */
                             nF,
                             resolution,
                             wingBreadth,
                             numLayers,
                             
                             out_d,
                             
                             T_d,
                             P_d,
                             
                             L[mol].iso,
                             L[mol].Vnn,
                             L[mol].Snn_ref,
                             L[mol].Yair,
                             L[mol].Yself,
                             L[mol].En,
                             L[mol].n,
                             L[mol].d,
                             &(N_d[(L[mol].mol-1) * atmosData->npfull]),
                             Z_d,                             
                             &(PS_d[(L[mol].mol-1) * atmosData->npfull]), /* offsets [mol] into idx@t,lat,lon,
                                                                           * so PS should be t, lat, lon, mol, z*/
                             &(LinesBuf_d[s]),
                             &(OpticsBuf_d[s]) );      
      }
    }    
  }

  /* block all streams */
  HANDLE_ERROR( cudaDeviceSynchronize() );

  /* free tmp buffers */
  for (s=0; s<*nStreams; ++s){
    device_optics_free(&(LinesBuf_d[s]),
                       &(OpticsBuf_d[s]) );
  }
  
  /* d2h */
  printf("Memcpy opticalDepth (out) to host..");
  HANDLE_ERROR( cudaMemcpy(out, out_d, numLayers*nF*sizeof(REAL_t),cudaMemcpyDeviceToHost) );
  printf("..done\n");
  
  device_atmos_free(T_d,
                    P_d,
                    N_d,
                    PS_d);
    
  HANDLE_ERROR(cudaFree(out_d));
  
  return EXIT_SUCCESS;
}

#else
int device_launch(int* nStreams,
                  void** streams,
                  const unsigned int numMols,
                  const char* const molFnames[],
                  const REAL_t loWn,
                  /* const REAL_t hiWn, */
                  const unsigned int nF,
                  const REAL_t resolution,
                  const unsigned int wingBreadth,
                  radiationOutputFields_t* atmosData,
                  REAL_t* const out,
                  int time,
                  int lat,
                  int lon)
{
  /* shut these wwarnings down */
  (void) nStreams;
  (void) streams;
  (void) numMols;
  (void) molFnames;
  (void) loWn;
  (void) nF;
  (void) resolution;
  (void) wingBreadth;
  (void) atmosData;
  (void) out;
  (void) time;
  (void) lat;
  (void) lon;
  printf("\n\nYou've not compiled with NVCC, device_launch does nothing...\n\n");
  return -1;
}
#endif
  
  
int main(int argc, char* argv[]){
  int world_size=-1;
  int world_rank=-1;
#ifdef MPI_ENABLED  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
#endif
  
  int time;
  unsigned int lat;
  unsigned int compute_lat_beg;
  unsigned int compute_lat_end;
  unsigned int lon;
  unsigned int compute_lon_beg;
  unsigned int compute_lon_end;
  unsigned int mol;
  size_t idx;
  REAL_t* out=NULL;
  
  /**************** parsing stuff ********************/

  /* where my args at */
  struct arguments arguments;
  /* defaults */
  arguments.silent = 0;
  arguments.verbose = 0;
  arguments.device=0;
  arguments.mpi=0;
  arguments.host=0;
  arguments.nhitfiles=0;
  arguments.nmolConc=0;
  arguments.nmolConcOver=0;
  arguments.atmos = NULL;
  static char default_output_fname[] ="didyouforgettospecifyoutfile.nc";
  arguments.output_file = default_output_fname;
  arguments.wingBreadth = 25;
  arguments.ctm = 0;
  /* model bounds */\
  arguments.w = 1;
  arguments.W = 3000;
  arguments.t = 0;
  arguments.T = 0;
  arguments.res = 1.0;
  /* molecular concentrations */
  arguments.h2o = 0;  /* default to nc after parsing */
  arguments.co2 = 0;
  arguments.o3 = 0;    /* default to nc after parsing */
  arguments.n2o = 0;
  arguments.co = 0;
  arguments.ch4 = 0;
  arguments.o2 = 0;  
  /* go */
  argp_parse( &argp,argc, argv, 0,0, &arguments);
  
  /* default default overrides */
  if(  arguments.h2o == 0 )
  {    
    arguments.h2o = -1;
    arguments.nmolConcOver++;
  }
  if(  arguments.o3 == 0 )
  {
    arguments.o3 = -1;
    arguments.nmolConcOver++;
  }

  /**************** ends parsing stuff ***********************/
    
  /**************** setup inputs and space for output **********************/

  const unsigned int nF = (arguments.W-arguments.w)/arguments.res + 1;
  
  if(arguments.device!=0 && arguments.host!=0)
  {
    fprintf(stderr, "Failure, Ambiguous launch type!  Specifiy singular flag --host/--device or neither flag to just default to device 0.\n");
    exit(EXIT_FAILURE);
  }
  else if(arguments.mpi!=0 && (arguments.device==0 && arguments.host==0) )
  {
    fprintf(stderr, "Failure, when specifying mpi you must specify number of devices per node or --host.\n");
    exit(EXIT_FAILURE);
  };

  const int launchType = arguments.host==1 ? 0 : 1 ;  /* launchType host is 0, device is 1, just a helper var */

#ifdef MPI_ENABLED
  const int device_number = world_size % arguments.device;
  assert( device_number>=0 );  /* if this ever trips try (a+n) % n */
  if( arguments.mpi!=0 && world_size <1)
  {
    fprintf(stderr,"You have invoked with --mpi but not used an mpirun style executer.  Confused\n");
    exit(EXIT_FAILURE);
  }
#else
  if(arguments.mpi!=0)
  {
    fprintf(stderr, "You have not built with MPI or set MPI_ENABLED\n");
    exit(EXIT_SUCCESS);
  }      
#endif
  
  if (launchType == 1)
  {
#ifdef __NVCC__
    const int device_number = arguments.device;
    HANDLE_ERROR( cudaSetDevice(device_number) );
#endif
  }
  
  char *atmosFile = arguments.atmos;
  
  radiationOutputFields_t atmosData;
  getAndSetAtmosFieldsFromFile(atmosFile, &atmosData);
  
  const unsigned int numLayers = atmosData.npfull;

  if(arguments.nmolConc!=arguments.nhitfiles)
  {
    fprintf(stderr,"Warning, Number of hitfiles (%d) does not match number of prescribed concentrations (%d), checking for sane overrides...\n",arguments.nhitfiles, arguments.nmolConc );
    if(arguments.nmolConc+arguments.nmolConcOver == arguments.nhitfiles)
    {
      fprintf(stderr,"\t...found %d overrides, okay computer.\n", arguments.nmolConcOver);
    }
    else
    {
      exit(EXIT_FAILURE);
    }
  }
  const unsigned int nMols = arguments.nhitfiles;
  char** hitFnameList = arguments.hitfiles;
  printf("\nSubmitted %d molecules.\n",nMols);

  /* init output file */
  int ncid;
  int varid;
  char *OUTPUT_FNAME=NULL;
  compute_lat_beg = 0;
  compute_lat_end = atmosData.nlat;
  compute_lon_beg = 0;
  compute_lon_end = atmosData.nlon;

  /* output file and compute setup is very different for mpi */
  if(arguments.mpi!=0)
  {
    OUTPUT_FNAME = (char*)malloc(strlen(arguments.output_file) + 9);
    if(OUTPUT_FNAME ==NULL)
    {
      fprintf(stderr, "Malloc failed for %zu bytes of OUTPUTFNAME!\n",
              strlen(arguments.output_file) + 9);
      exit(EXIT_FAILURE);
    }
    
    lat = atmosData.nlat / world_size;
    if(lat*world_size != atmosData.nlat)
    {
      fprintf(stderr, "Warning:"
	      "\n\tSpecified %d global lats across ranks=%zu yields between %d and %d lats per rank.  "
	      "\n\tThis will result in idle hardware, suggest a different work share.\n",
	      world_size,
	      atmosData.nlat,
	      lat,
	      lat+1);
    }
    compute_lat_beg = world_rank*lat;
    compute_lat_end = compute_lat_beg+lat;
    if( compute_lat_end>atmosData.nlat)
    {
      compute_lat_end = atmosData.nlat;
    }
    compute_lon_beg = 0;
    compute_lon_end = atmosData.nlon;
    /* copy existing name
     * cat .rankN */
    sprintf(OUTPUT_FNAME, "%s.rank%d", arguments.output_file, world_rank);
  }
  else
  {
    OUTPUT_FNAME = arguments.output_file;
  }
  fprintf(stderr,"Opening output file %s.\n",OUTPUT_FNAME);
  openOpticalDepthOutput(&ncid, &varid, OUTPUT_FNAME,
                         compute_lat_end - compute_lat_beg,
                         compute_lon_end - compute_lon_beg,
                         numLayers,
                         nF);
#ifdef __NVCC__
  int nstreams=-1;
  cudaStream_t* streams=NULL;
#endif

  /* setup hitran lines */
  RefLinePtrs_t HitLines[nMols];
  
  /* get fname and parse in lines */
  RefLine_flags_t flags= {((unsigned int) -1),1,0}; /* host cuda malloc default, host=True, device=false */
  time = arguments.T - arguments.t + 1;
  for(mol=0; mol<nMols; ++mol){
      HitLines[mol] = parseHITRANfile(hitFnameList[mol],flags,arguments.w,arguments.W);
      /* checkMolConfig(&arguments, HitLines[mol].mol, atmosData.PS, time ,atmosData.nlat, atmosData.nlon, numLayers ); */
      checkMolConfig(&arguments, HitLines[mol].mol, &atmosData, time);
  }

  for( time=arguments.t ; time<=arguments.T ; ++time)
  {
    for(lat=compute_lat_beg; lat<compute_lat_end; ++lat)
    {
      for( lon=compute_lon_beg; lon<compute_lon_end; ++lon)
      {
        if ( launchType == 0){
          if(out == NULL)  /* malloc if needed, otherwise pass */
          {
            out = (REAL_t*)calloc(nF*numLayers,sizeof(REAL_t));
          }
          
          host_launch( nMols ,
                       /* hitFnameList, */
                       HitLines,
                       ((REAL_t)arguments.w),
                       /* ((REAL_t)arguments.W), */
                       nF,
                       arguments.res,
                       arguments.wingBreadth,
                       &atmosData,
                       out,
                       time,
                       lat,
                       lon);  
        }
        else if (launchType==1){
          
#ifdef __NVCC__
          if(out == NULL)  /* malloc if needed, otherwise pass */
          {
            HANDLE_ERROR( cudaHostAlloc(&out , nF*numLayers*sizeof(REAL_t) , cudaHostAllocDefault ) );
          }
          
          device_launch( &nstreams,
                         &streams,
                         nMols,
                         /* hitFnameList, */
                         HitLines,
                         ((REAL_t)arguments.w),
                         /* ((REAL_t)arguments.W), */
                         nF,
                         arguments.res,
                         arguments.wingBreadth,
                         &atmosData,
                         out,
                         time,
                         lat,
                         lon);
#else
          fprintf(stderr,"\n Requested cuda launch type (%d), but compiled host only... Aborting.\n",launchType);
          exit(1);
#endif    
        }
        else {
          fprintf(stderr,"\n Unkown launch type (%d) requested. Aborting.\n",launchType);
          exit(1);
        }

        if(arguments.ctm==1)
        {
          printf("Computing Continuum\n");
          assert(arguments.res==1.);  /* presently the continuum code is only safe for widths of one wavenumber */
          idx = (time*atmosData.nlon*atmosData.nlat + lat*atmosData.nlon + lon) * atmosData.npfull;
          /* dbg print */
          printf("%zu: T=%g P=%g DELTAZ=%g PS[%zu]=%g \n",
                 idx,
                 atmosData.T[idx],
                 atmosData.P[idx],
                 atmosData.DELTAZ[idx],
                 NUM_MOL*idx + H2O*atmosData.npfull,
                 atmosData.PS[NUM_MOL*idx + H2O*atmosData.npfull] );
          /* getchar(); */
          get_CTM(out,
                  &(atmosData.T[idx]),
                  &(atmosData.P[idx]),
                  &(atmosData.DELTAZ[idx]),
                  &(atmosData.PS[NUM_MOL*idx + H2O*atmosData.npfull ]),
                  nF,
                  atmosData.npfull);
        }

        fprintf(stderr, "Writing hyperslab of %d samples "
                "@{t=%d, lat=%d, lon=%d, layers=0:%d} to output file %s \n",
              nF, time, lat, lon, numLayers, OUTPUT_FNAME);
        writeOpticalDepthOutputByColumn(ncid,
					varid,
					time,
					lat-compute_lat_beg,
					lon-compute_lon_beg,
					numLayers,
					nF,
					out);

        memset(out, 0, nF*numLayers*sizeof(REAL_t));
        
      }  /* nlat */
    }  /* nlon */
  }  /* time */
  
  fprintf(stderr, "Closing output file %s\n", OUTPUT_FNAME);
  closeOpticalDepthOutput(ncid);
  
  /* cleanup all, this is dirty,  into earlier stage later */
  for(mol=0; mol<nMols; ++mol){
    freeHost(HitLines[mol],flags);
  }
  
  
///wrap up something like this in a function, then key off of outputfile extension for csv output
/* #undef WRITEOUT */
/*   ///#define WRITEOUT */
/* #ifdef WRITEOUT */
/*   printf("\n\n\t Attempting Result Write.\n\n"); */
/*   /\* output *\/ */
/*   REAL_t wv; */
/*   FILE* ofp; */
/*   ofp = fopen("opticaldepth.testout.csv","w"); */
/*   if(ofp==NULL){ */
/*     fprintf(stderr,"\nopening output file for writing failed, aborting.\n"); */
/*     exit(1); */
/*   } */
/*   /\* header *\/ */
/*   fprintf(ofp, "z,wavenumber,val\n"); */
/*   /\* data *\/ */
/*   for(unsigned int l=0; l<numLayers; ++l){ */
/*     /\* for(unsigned int iter=0; iter<nF; iter+=(((double)1)/arguments.res) ){ /\\* output every one wavenumber *\\/ *\/ */
/*     for(unsigned int iter=0; iter<nF; ++iter ){ /\* output every sample *\/ */
/*       wv = iter*arguments.res + arguments.w; */
/*       assert(wv <= arguments.W); */
/*       fprintf(ofp,"%u %.17f %.17f\n", l, wv , out[ l*nF + iter ] ); */
/*     } */
/*   } */
  
/*   if( fclose(ofp)!=0 ){ */
/*     fprintf(stderr,"\nclosing output file for writing failed, aborting.\n"); */
/*     exit(1); */
/*   } */
/* #endif */

  /* cleanup */
  if (launchType==1){
#ifdef __NVCC__
    cudaFreeHost(out);
#endif
  }
  else{
    free(out);
  }

#ifdef MPI_ENABLED
  MPI_Finalize();
#endif
  
  return EXIT_SUCCESS;
}

