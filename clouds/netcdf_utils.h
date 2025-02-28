#ifndef NETCDF_UTILS_H
#define NETCDF_UTILS_H

#include "netcdf.h"


/* @brief Crashes if any netCDF errors are detected.*/
void netcdf_catch(int const err);


/* @brief Opens netCDF dataset.*/
int open_dataset(char const * path);


/* @brief Closes netCDF dataset.*/
void close_dataset(int const ncid);


/* @brief Reads an attribute from netCDF dataset.*/
void read_attribute(int const ncid, char const * variable, char const * name,
                    nc_type const type, void * buffer);


/* @brief Read the length of a dimension from a netCDF dataset.*/
void read_dimlen(int const ncid, char const * name, int * length);


/* @brief Reads variable from netCDF dataset.*/
void read_variable(int const ncid, char const * name, void ** buffer,
                   nc_type const type, size_t const * start, size_t const * count);


#endif
