#ifndef OPTICS_UTILS_H
#define OPTICS_UTILS_H


/* @brief Optical properties.*/
typedef struct OpticalProperties
{
    double * bands; /*Band center wavenumbers [cm-1] (band).*/
    double * band_limits; /*Band edge wavenumbers [cm-1] (band, 2).*/
    double * extinction_coefficient; /*Extinction coefficient (band).*/
    int num_bands; /*Number of bands.*/
    double * single_scatter_albedo; /*Single-scatter albedo (band).*/
    double * asymmetry_factor; /*Asymmetry factor (band).*/
} OpticalProperties_t;


/* @brief Constructs an OpticalProperties object.*/
void construct_optics(OpticalProperties_t * self, int const num_bands, double const * bands,
                      double const * band_limits);


/* @brief Destructs an OpticalProperties object.*/
void destruct_optics(OpticalProperties_t * self);


/* @brief Creates a new OpticalProperties object on a new set of bands from the current
          OpticalProperties object.*/
void thick_average(OpticalProperties_t const self, OpticalProperties_t * optics,
                   int const * starting_band, int const * ending_band);


#endif
