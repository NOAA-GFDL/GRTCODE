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
#include <assert.h>
#include <netcdf.h>
#include "outputNetcdfSpec.h"


/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define NCERR(e) {fprintf(stderr, "Error: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}


void closeOpticalDepthOutput(int ncid){
  int retval;
  /* Close the file. This frees up any internal netCDF resources
   * associated with the file, and flushes any buffers. */
  if ((retval = nc_close(ncid)))
    NCERR(retval);
}

void openOpticalDepthOutput(int* ncid,
                            int* varid,
                            char FNAME[],
                            const size_t nlat,
                            const size_t nlon,
                            const size_t nlayers,
                            const size_t nF )
{
  int retval;
  int t_dimid;
  int lat_dimid;
  int lon_dimid;
  int lay_dimid;
  int f_dimid;
  const int ndims=5;
  int dimids[ndims];
      
  /* Create the file.
   * The NC_CLOBBER parameter tells netCDF to overwrite this file, if it already exists.  */
  if ((retval = nc_create(FNAME, NC_CLOBBER, ncid)))
    NCERR(retval);

  /* NOFILL does not prefill file, so avoids extraneous writing */
  ncsetfill(*ncid, NC_NOFILL); 

  /* Define the dimensions. NetCDF will hand back an ID for each. */
  if ((retval = nc_def_dim(*ncid, "time", NC_UNLIMITED, &t_dimid)))
    NCERR(retval);
  if ((retval = nc_def_dim(*ncid, "lat", nlat, &lat_dimid)))
    NCERR(retval);
  if ((retval = nc_def_dim(*ncid, "lon", nlon, &lon_dimid)))
    NCERR(retval);
  if ((retval = nc_def_dim(*ncid, "pfull", nlayers, &lay_dimid)))
    NCERR(retval);
  if ((retval = nc_def_dim(*ncid, "wavenumber", nF, &f_dimid)))
    NCERR(retval);

  dimids[0] = t_dimid;
  dimids[1] = lat_dimid;
  dimids[2] = lon_dimid;
  dimids[3] = lay_dimid;
  dimids[4] = f_dimid;

  /* Define the variable. */
  if ((retval = nc_def_var(*ncid, "OpticalDepth", NC_FLOAT, ndims,
                           dimids, varid)))
    NCERR(retval);
  
  /* End define mode. This tells netCDF we are done defining
   * metadata. */
  if ((retval = nc_enddef(*ncid)))
    NCERR(retval)
}
                             

void writeOpticalDepthOutputByColumn(const int ncid,
                                    const int varid,
                                    const int t,                                    
                                    const int lat,
                                    const int lon,
                                    const int nlayers,
                                    const int nF,
                                    float* spectra)
{
  int retval;
  const int ndims = 5;
  size_t start[ndims];
  size_t count[ndims];

  /* a column in time */
  count[0] = 1;  /* 1 time */
  count[1] = 1;  /* 1 lat */
  count[2] = 1;  /* 1 lon */
  /* is composed of  */
  count[3] = nlayers;  /* layers in the column */
  count[4] = nF;       /* samples per layer */

  /* the column we inted to write is at */
  start[0] = t;
  start[1] = lat;
  start[2] = lon;
  start[3] = 0;  /* zeroth layer */
  start[4] = 0;  /* zeroth sample */

  if ((retval = nc_put_vara_float(ncid, varid, start, count,
                                  spectra)))
    NCERR(retval);

}
      
