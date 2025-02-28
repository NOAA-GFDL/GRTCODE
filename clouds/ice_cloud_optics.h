#ifndef ICE_CLOUD_OPTICS_H
#define ICE_CLOUD_OPTICS_H

#include "optics_utils.h"


/* @brief Ice water cloud optics parameterizations.*/
typedef struct IceCloudOptics
{
    double * a; /*a parameters from equations 4a/5a (band, 6).*/
    double * b; /*b parameters from equations 4b/5b (band, 6).*/
    double * band_limits; /*Parameterization band limits [cm-1] (bands, 2).*/
    double * bands; /*Parameterization band centers [cm-1] (band).*/
    double * c; /*c parameters from equations 4c/5c (radius, band, 6).*/
    int last_ir_band; /*Index of last infrared (longwave) band.*/
    int num_bands; /*Number or bands.*/
    int num_order; /*Number of polynomial orders.*/
    int num_radius_bins; /*Number of radius bins.*/
    double * radii; /*Radius bins [micron] for the parameterization (radius, 2).*/
} IceCloudOptics_t;


/* @brief Constructs a IceCloudOptics object.*/
void construct_ice_optics(IceCloudOptics_t * self, char const * path);


/* @brief Constructs a IceCloudOptics object.*/
void destruct_ice_optics(IceCloudOptics_t * self);


/* @brief Calculates cloud optics.*/
void calculate_ice_optics(IceCloudOptics_t const self, double const ice_concentration,
                          double const equivalent_radius, double const scale_factor,
                          double const temperature, OpticalProperties_t * optical_properties);


#endif
