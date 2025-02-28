/* @brief NetCDF utilities.*/
#include <stdio.h>
#include <stdlib.h>
#include "netcdf_utils.h"
#include "netcdf.h"


/* @brief Crashes if any netCDF errors are detected.*/
void netcdf_catch(int const err)
{
    if (err != NC_NOERR)
    {
        fprintf(stderr, "%s", nc_strerror(err));
        exit(1);
    }
    return;
}


/* @brief Opens netCDF dataset.*/
int open_dataset(char const * path)
{
    int ncid;
    netcdf_catch(nc_open(path, NC_NOWRITE, &ncid));
    return ncid;
}


/* @brief Closes netCDF dataset.*/
void close_dataset(int const ncid)
{
    netcdf_catch(nc_close(ncid));
    return;
}


/* @brief Reads an attribute from netCDF dataset.*/
void read_attribute(int const ncid, char const * variable, char const * name,
                    nc_type const type, void * buffer)
{
    int varid;
    netcdf_catch(nc_inq_varid(ncid, variable, &varid));
    if (type == NC_INT)
    {
        netcdf_catch(nc_get_att_int(ncid, varid, name, (int *)buffer));
    }
    else if (type == NC_FLOAT)
    {
        netcdf_catch(nc_get_att_float(ncid, varid, name, (float *)buffer));
    }
    else if (type == NC_DOUBLE)
    {
        netcdf_catch(nc_get_att_double(ncid, varid, name, (double *)buffer));
    }
    else if (type == NC_CHAR)
    {
        netcdf_catch(nc_get_att_text(ncid, varid, name, (char *)buffer));
    }
    else
    {
        fprintf(stderr, "non-supported netcdf type.");
        exit(1);
    }
    return;
}


/* @brief Read the length of a dimension from a netCDF dataset.*/
void read_dimlen(int const ncid, char const * name, int * length)
{
    int dimid;
    netcdf_catch(nc_inq_dimid(ncid, name, &dimid));
    size_t size;
    netcdf_catch(nc_inq_dimlen(ncid, dimid, &size));
    *length = (int)size;
    return;
}


/* @brief Reads variable from netCDF dataset.*/
void read_variable(int const ncid, char const * name, void ** buffer,
                   nc_type const type, size_t const * start, size_t const * count)
{
    int varid;
    netcdf_catch(nc_inq_varid(ncid, name, &varid));
    int ndims;
    netcdf_catch(nc_inq_varndims(ncid, varid, &ndims));
    size_t corner[ndims];
    if (start == NULL)
    {
        int i;
        for (i=0; i<ndims; ++i)
        {
             corner[i] = 0;
        }
    }
    else
    {
        int i;
        for (i=0; i<ndims; ++i)
        {
            corner[i] = start[i];
        }
    }
    size_t sizes[ndims];
    if (count == NULL)
    {
        int dimids[ndims];
        netcdf_catch(nc_inq_vardimid(ncid, varid, dimids));
        int i;
        for (i=0; i<ndims; ++i)
        {
            netcdf_catch(nc_inq_dimlen(ncid, dimids[i], &(sizes[i])));
            sizes[i] -= corner[i];
        }
    }
    else
    {
        int i;
        for (i=0; i<ndims; ++i)
        {
            sizes[i] = count[i];
        }
    }
    size_t total_size = 1;
    int i;
    for (i=0; i<ndims; ++i)
    {
        total_size *= sizes[i];
    }
    size_t num_bytes;
    if (type == NC_INT)
    {
        num_bytes = sizeof(int);
    }
    else if (type == NC_FLOAT)
    {
        num_bytes = sizeof(float);
    }
    else if (type == NC_DOUBLE)
    {
        num_bytes = sizeof(double);
    }
    else
    {
        fprintf(stderr, "non-supported netcdf type.");
        exit(1);
    }
    void * data = malloc(num_bytes*total_size);
    if (type == NC_INT)
    {
        netcdf_catch(nc_get_vara_int(ncid, varid, corner, sizes, (int *)data));
    }
    else if (type == NC_FLOAT)
    {
        netcdf_catch(nc_get_vara_float(ncid, varid, corner, sizes, (float *)data));
    }
    else if (type == NC_DOUBLE)
    {
        netcdf_catch(nc_get_vara_double(ncid, varid, corner, sizes, (double *)data));
    }
    *buffer = data;
    return;
}
