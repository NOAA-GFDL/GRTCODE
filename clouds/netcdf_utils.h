#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H


#include "device.h"
#include "netcdf.h"


/* @brief Opens netCDF dataset.*/
int open_dataset(char const * path, int * ncid);


/* @brief Closes netCDF dataset.*/
int close_dataset(int const ncid);


/* @brief Reads an attribute from netCDF dataset.*/
int read_attribute(int const ncid, char const * variable, char const * name,
                   nc_type const type, void * buffer);


/* @brief Read the length of a dimension from a netCDF dataset.*/
int read_dimlen(int const ncid, char const * name, int * length);


/* @brief Reads variable from netCDF dataset.*/
int read_variable(int const ncid, char const * name, void ** buffer,
                  nc_type const type, Device_t const device,
                  size_t const * start, size_t const * count);


#endif
