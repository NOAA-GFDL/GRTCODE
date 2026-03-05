/* @brief Incomplete beta distribution utilties.*/
#include <stdlib.h>
#include "incomplete_beta.h"
#include "netcdf_utils.h"


/* @brief Constructs an IncompleteBeta object.*/
void construct_beta(IncompleteBeta_t * self, char const * path)
{
    int ncid = open_dataset(path);
    read_dimlen(ncid, "p", &(self->num_shape));
    read_dimlen(ncid, "x", &(self->num_x));
    read_variable(ncid, "p", (void **)&(self->p), NC_INT, NULL, NULL);
    read_variable(ncid, "q", (void **)&(self->q), NC_INT, NULL, NULL);
    read_variable(ncid, "x", (void **)&(self->x), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "data", (void **)&(self->y), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "inverse", (void **)&(self->y_inverse), NC_DOUBLE, NULL, NULL);
    close_dataset(ncid);
    return;
}


/* @brief Destructs an IncompleteBeta object.*/
void destruct_beta(IncompleteBeta_t * self)
{
    free(self->p);
    free(self->q);
    free(self->x);
    free(self->y);
    free(self->y_inverse);
    return;
}


/* @brief One-dimensional linear interpolation.*/
static double interp(int const num_x, double const * x, double const * y, double const newx)
{
    int i;
    for (i=1; i<(num_x - 1); ++i)
    {
        if (x[i] > newx)
        {
            break;
        }
    }
    double const m = (y[i] - y[i-1])/(x[i] - x[i-1]);
    double const b = y[i] - m*x[i];
    return m*newx + b;
}


/* @brief Calculates the inverse of the incomplete beta function.*/
double beta_inverse(IncompleteBeta_t const self, int const p, int const q, double const x)
{
    int const offset = ((q - 1)*self.num_shape + (p - 1))*self.num_x;
    return interp(self.num_x, self.x, &(self.y_inverse[offset]), x);
}


/* @brief Calculates incomplete beta function.*/
double beta_value(IncompleteBeta_t const self, int const p, int const q, double const x)
{
    int const offset = ((q - 1)*self.num_shape + (p - 1))*self.num_x;
    return interp(self.num_x, self.x, &(self.y[offset]), x);
}