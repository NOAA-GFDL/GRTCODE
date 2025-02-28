#include <stdlib.h>
#include <string.h>
#include "optics_utils.h"


/* @brief Constructs an OpticalProperties object.*/
void construct_optics(OpticalProperties_t * self, int const num_bands, double const * bands,
                      double const * band_limits)
{
    self->num_bands = num_bands;
    self->bands = malloc(sizeof(double)*num_bands);
    memcpy(self->bands, bands, sizeof(double)*num_bands);
    self->band_limits = malloc(sizeof(double)*num_bands*2);
    memcpy(self->band_limits, band_limits, sizeof(double)*num_bands*2);
    self->extinction_coefficient = malloc(sizeof(double)*num_bands);
    self->single_scatter_albedo = malloc(sizeof(double)*num_bands);
    self->asymmetry_factor = malloc(sizeof(double)*num_bands);
    return;
}


/* @brief Destructs an OpticalProperties object.*/
void destruct_optics(OpticalProperties_t * self)
{
    free(self->bands);
    free(self->band_limits);
    free(self->extinction_coefficient);
    free(self->single_scatter_albedo);
    free(self->asymmetry_factor);
    return;
}


/* @brief Linearly interpolates.*/
static double interp(double const * x, double const * y, double const newx)
{
    double const m = (y[1] - y[0])/(x[1] - x[0]);
    double const b = y[0] - m*x[0];
    return m*newx + b;
}


/* @brief Creates a new OpticalProperties object on a new set of bands from the current
          OpticalProperties object.*/
void thick_average(OpticalProperties_t const self, OpticalProperties_t * optics,
                   int const * starting_band, int const * ending_band)
{
    int const a = starting_band == NULL ? 0 : *starting_band - 1;
    int const b = ending_band == NULL ? self.num_bands - 1 : *ending_band - 1;

    /*Linearly interpolate for now.*/
    int i;
    for (i=0; i<optics->num_bands; ++i)
    {
        if (self.band_limits[a*2] < optics->bands[i]) break;
    }
    if (i > 0)
    {
        /*Use value from the first band in all new bands that are less than the first band's
          lower limit.*/
        int j;
        for (j=0; j<i; ++j)
        {
            optics->extinction_coefficient[j] = self.extinction_coefficient[a];
            optics->single_scatter_albedo[j] = self.single_scatter_albedo[a];
            optics->asymmetry_factor[j] = self.asymmetry_factor[a];
        }
    }
    int c = a;
    int j;
    for(j=i; j<optics->num_bands; ++j)
    {
        if (optics->bands[j] > self.band_limits[2*b + 1]) break;
        int k;
        for (k=c; k<=b; ++k)
        {
            if (self.bands[k] > optics->bands[j]) break;
        }
        c = k;
        if (k == a)
        {
            if (optics->bands[j] <= self.bands[k])
            {
                /*Inside the first band, but left of the band center.*/
                optics->extinction_coefficient[j] = self.extinction_coefficient[k];
                optics->single_scatter_albedo[j] = self.single_scatter_albedo[k];
                optics->asymmetry_factor[j] = self.asymmetry_factor[k];
                continue;
            }
            k = a + 1;
        }
        if (k == b + 1)
        {
            k = b;
            if (optics->bands[j] >= self.bands[k])
            {
                /*Inside the last band, but right of the band center.*/
                optics->extinction_coefficient[j] = self.extinction_coefficient[k];
                optics->single_scatter_albedo[j] = self.single_scatter_albedo[k];
                optics->asymmetry_factor[j] = self.asymmetry_factor[k];
                continue;
            }
        }
        optics->extinction_coefficient[j] = interp(&(self.bands[k-1]),
                                                   &(self.extinction_coefficient[k-1]),
                                                   optics->bands[j]);
        optics->single_scatter_albedo[j] = interp(&(self.bands[k-1]),
                                                  &(self.single_scatter_albedo[k-1]),
                                                  optics->bands[j]);
        optics->asymmetry_factor[j] = interp(&(self.bands[k-1]),
                                             &(self.asymmetry_factor[k-1]),
                                             optics->bands[j]);
    }
    if (j <= optics->num_bands - 1)
    {
        /*Use value from the last band in all new bands that are greater than the last
          band's upper limit.*/
        int k;
        for (k=j; k<optics->num_bands; ++k)
        {
            optics->extinction_coefficient[k] = self.extinction_coefficient[b];
            optics->single_scatter_albedo[k] = self.single_scatter_albedo[b];
            optics->asymmetry_factor[k] = self.asymmetry_factor[b];
        }
    }
    return;
}
