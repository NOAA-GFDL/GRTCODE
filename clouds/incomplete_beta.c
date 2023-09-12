/* @brief Incomplete beta distribution utilties.*/
#include <stdlib.h>
#include "debug.h"
#include "incomplete_beta.h"
#include "netcdf.h"
#include "netcdf_utils.h"


#include <stdio.h>


/* @brief Constructs an IncompleteBeta object.*/
EXTERN int create_beta(IncompleteBeta_t * self, char const * path,
                       Device_t const * const device)
{
    not_null(self);
    not_null(path);
    not_null(device);
    self->device = *device;
    int ncid;
    catch(open_dataset(path, &ncid));
    catch(read_dimlen(ncid, "p", &(self->num_shape)));
    catch(read_dimlen(ncid, "x", &(self->num_x)));
    nc_type var_type;
    if (sizeof(fp_t) == sizeof(double))
    {
        var_type = NC_DOUBLE;
    }
    else
    {
        var_type = NC_FLOAT;
    }
    catch(read_variable(ncid, "p", (void **)&(self->p), NC_INT, self->device, NULL, NULL));
    catch(read_variable(ncid, "q", (void **)&(self->q), NC_INT, self->device, NULL, NULL));
    catch(read_variable(ncid, "x", (void **)&(self->x), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "data", (void **)&(self->y), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "inverse", (void **)&(self->y_inverse), var_type, self->device, NULL, NULL));
    catch(close_dataset(ncid));
    return GRTCODE_SUCCESS;
}


/* @brief Destructs an IncompleteBeta object.*/
EXTERN int destroy_beta(IncompleteBeta_t * self)
{
    gfree(self->p, self->device);
    gfree(self->q, self->device);
    gfree(self->x, self->device);
    gfree(self->y, self->device);
    gfree(self->y_inverse, self->device);
    return GRTCODE_SUCCESS;
}


/* @brief One-dimensional linear interpolation.*/
HOST DEVICE static fp_t interp(int const num_x, fp_t const * x,
                               fp_t const * y, fp_t const newx)
{
    int i;
    for (i=1; i<(num_x - 1); ++i)
    {
        if (x[i] > newx)
        {
            break;
        }
    }
    fp_t const m = (y[i] - y[i-1])/(x[i] - x[i-1]);
    fp_t const b = y[i] - m*x[i];
    return m*newx + b;
}


/* @brief Calculates the inverse of the incomplete beta function.*/
HOST DEVICE fp_t beta_inverse(int const num_shape, int const num_x, 
                              fp_t const * x, fp_t const * y_inverse, int const p,
                              int const q, fp_t const z)
{
    int const offset = ((q - 1)*num_shape + (p - 1))*num_x;
    return interp(num_x, x, &(y_inverse[offset]), z);
}


/* @brief Calculates incomplete beta function.*/
HOST DEVICE fp_t beta_value(int const num_shape, int const num_x,
                            fp_t const * x, fp_t const * y, int const p,
                            int const q, fp_t const z)
{
    int const offset = ((q - 1)*num_shape + (p - 1))*num_x;
    return interp(num_x, x, &(y[offset]), z);
}
