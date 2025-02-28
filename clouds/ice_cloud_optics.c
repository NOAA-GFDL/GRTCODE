/* @brief Ice water cloud optics parameterizations from doi: 10.1175/2009JCLI2844.1*/
#include <stdlib.h>
#include "ice_cloud_optics.h"
#include "netcdf_utils.h"
#include "optics_utils.h"


/* @brief Donner parameterization.*/
static double ice_particle_size(double const temperature)
{
    double const tfreeze = 273.16;
    if (temperature > tfreeze - 25.)
    {
        return 100.6;
    }
    else if (temperature > tfreeze - 30.)
    {
        return 80.8;
    }
    else if (temperature > tfreeze - 35.)
    {
        return 93.5;
    }
    else if (temperature > tfreeze - 40.)
    {
        return 63.9;
    }
    else if (temperature > tfreeze - 45.)
    {
        return 42.5;
    }
    else if (temperature > tfreeze - 50.)
    {
        return 39.9;
    }
    else if (temperature > tfreeze - 55.)
    {
        return 21.6;
    }
    else
    {
        return 20.2;
    }
}


/* @brief Calculates cloud optics.*/
static void optics(IceCloudOptics_t const self, double const ice_concentration,
                   double const equivalent_radius, double const scale_factor,
                   double const temperature, int const band,
                   double * extinction_coefficient, double * single_scatter_albedo,
                   double * asymmetry_factor)
{
    double radius = equivalent_radius > 0. ? equivalent_radius : ice_particle_size(temperature);
    int r;
    for (r=1; r<self.num_radius_bins; ++r)
    {
        if (self.radii[2*r] > radius) break;
    }
    r -= 1;
    double const min_radius = 13.;
    radius = min_radius > scale_factor*radius ? min_radius : scale_factor*radius;
    double d[self.num_order];
    d[0] = 1.;
    int i;
    for (i=1; i<self.num_order; ++i)
    {
        d[i] = d[i - 1]*radius;
    }
    double d_inv[self.num_order];
    for (i=0; i<self.num_order; ++i)
    {
        d_inv[i] = 1./d[i];
    }

    double asum = 0.;
    for (i=0; i<self.num_order; ++i)
    {
        asum += self.a[band*self.num_order + i]*d_inv[i];
    }
    *extinction_coefficient = ice_concentration*asum;
    if (band < self.last_ir_band)
    {
        double bsum = 0.;
        for (i=0; i<self.num_order; ++i)
        {
            bsum += self.b[band*self.num_order + i]*d_inv[i];
        }
        *single_scatter_albedo = 1. - (ice_concentration*bsum/(*extinction_coefficient));
    }
    else
    {
        double bsum = 0.;
        for (i=0; i<self.num_order; ++i)
        {
            bsum += self.b[band*self.num_order + i]*d[i];
        }
        *single_scatter_albedo = 1. - bsum;
    }
    double csum = 0.;
    for (i=0; i<self.num_order; ++i)
    {
        csum += self.c[(r*self.num_bands + band)*self.num_order + i]*d[i];
    }
    *asymmetry_factor = csum;
    return;
}


/* @brief Constructs a IceCloudOptics object.*/
void construct_ice_optics(IceCloudOptics_t * self, char const * path)
{
    int ncid = open_dataset(path);
    read_variable(ncid, "radius_bnds", (void **)&(self->radii), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "band_bnds", (void **)&(self->band_limits), NC_DOUBLE, NULL, NULL);
    read_dimlen(ncid, "band", &(self->num_bands));
    read_dimlen(ncid, "order5", &(self->num_order));
    read_dimlen(ncid, "radius", &(self->num_radius_bins));
    self->bands = malloc(sizeof(double)*self->num_bands);
    int i;
    for (i=0; i<self->num_bands; ++i)
    {
        self->bands[i] = 0.5*(self->band_limits[2*i] + self->band_limits[2*i + 1]);
    }
    read_attribute(ncid, "band_bnds", "last_IR_band", NC_INT, &(self->last_ir_band));
    read_variable(ncid, "a", (void **)&(self->a), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "b", (void **)&(self->b), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "c", (void **)&(self->c), NC_DOUBLE, NULL, NULL);
    close_dataset(ncid);
    return;
}


/* @brief Constructs a IceCloudOptics object.*/
void destruct_ice_optics(IceCloudOptics_t * self)
{
    free(self->radii);
    free(self->band_limits);
    free(self->bands);
    free(self->a);
    free(self->b);
    free(self->c);
    return;
}


/* @brief Calculates cloud optics.*/
void calculate_ice_optics(IceCloudOptics_t const self, double const ice_concentration,
                          double const equivalent_radius, double const scale_factor,
                          double const temperature, OpticalProperties_t * optical_properties)
{
    int i;
    for (i=0; i<self.num_bands; ++i)
    {
        optics(self, ice_concentration, equivalent_radius, scale_factor, temperature, i,
               &(optical_properties->extinction_coefficient[i]),
               &(optical_properties->single_scatter_albedo[i]),
               &(optical_properties->asymmetry_factor[i]));
    }
    return;
}
