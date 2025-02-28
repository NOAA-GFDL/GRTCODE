/* @brief Hu and Stamnes cloud liquid water parameterization from
          doi: 10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "hu_stamnes.h"
#include "netcdf_utils.h"
#include "optics_utils.h"


/* @brief Calculates cloud optics.*/
static void optics(HuStamnes_t const self, double const water_concentration,
                   double const equivalent_radius, int const band,
                   double * extinction_coefficient, double * single_scatter_albedo,
                   double * asymmetry_factor)
{
    double r = self.min_radius > equivalent_radius ? self.min_radius : equivalent_radius;
    r = self.max_radius < r ? self.max_radius : r;
    int i;
    for (i=1; i<self.num_radius_bins; ++i)
    {
        if (self.radii[i] > r) break;
    }
    i = (i - 1)*self.num_bands + band;
    double const m_to_km = 1.e-3;
    *extinction_coefficient = water_concentration*m_to_km*(self.a1[i]*
                              pow(r, self.b1[i]) + self.c1[i]); /*Equation 13.*/
    *single_scatter_albedo = 1. - (self.a2[i]*pow(r, self.b2[i]) + self.c2[i]); /*Equation 14.*/
    *asymmetry_factor = self.a3[i]*pow(r, self.b3[i]) + self.c3[i]; /*Equation 15.*/
    return;
}


/* @brief Constructs a HuStamnes object.*/
void construct_liquid_optics(HuStamnes_t * self, char const * path)
{
    int ncid = open_dataset(path);
    char name[128];
    memset(name, '\0', 128);
    read_attribute(ncid, "radius", "bounds", NC_CHAR, name);
    float radii[2];
    read_attribute(ncid, "radius", "valid_range", NC_FLOAT, radii);
    self->min_radius = (double)(radii[0]);
    self->max_radius = (double)(radii[1]);
    double * radius_bounds;
    read_variable(ncid, name, (void **)&radius_bounds, NC_DOUBLE, NULL, NULL);
    read_dimlen(ncid, "band", &self->num_bands);
    read_dimlen(ncid, "radius", &self->num_radius_bins);
    self->radii = malloc(sizeof(double)*(self->num_radius_bins + 1));
    int i;
    for (i=0; i<self->num_radius_bins; ++i)
    {
        self->radii[i] = radius_bounds[2*i];
    }
    self->radii[self->num_radius_bins] = radius_bounds[2*(self->num_radius_bins - 1) + 1];
    free(radius_bounds);
    read_variable(ncid, "band_bnds", (void **)&(self->band_limits), NC_DOUBLE,
                  NULL, NULL);
    read_attribute(ncid, "band_bnds", "last_IR_band", NC_INT, &self->last_ir_band);
    read_variable(ncid, "band", (void **)&(self->bands), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "a1", (void **)&(self->a1), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "a2", (void **)&(self->a2), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "a3", (void **)&(self->a3), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "b1", (void **)&(self->b1), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "b2", (void **)&(self->b2), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "b3", (void **)&(self->b3), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "c1", (void **)&(self->c1), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "c2", (void **)&(self->c2), NC_DOUBLE, NULL, NULL);
    read_variable(ncid, "c3", (void **)&(self->c3), NC_DOUBLE, NULL, NULL);
    close_dataset(ncid);
    return;
}


/* @brief Destructs a HuStamnes object.*/
void destruct_liquid_optics(HuStamnes_t * self)
{
    free(self->a1);
    free(self->a2);
    free(self->a3);
    free(self->band_limits);
    free(self->bands);
    free(self->b1);
    free(self->b2);
    free(self->b3);
    free(self->c1);
    free(self->c2);
    free(self->c3);
    free(self->radii);
    return;
}


/* @brief Calculates cloud optics.*/
void calculate_liquid_optics(HuStamnes_t const self, double const water_concentration,
                             double const equivalent_radius,
                             OpticalProperties_t * optical_properties)
{
    int i;
    for (i=0; i<self.num_bands; ++i)
    {
        optics(self, water_concentration, equivalent_radius, i,
               &(optical_properties->extinction_coefficient[i]),
               &(optical_properties->single_scatter_albedo[i]),
               &(optical_properties->asymmetry_factor[i]));
    }
    return;
}
