#ifndef HU_STAMNES_H
#define HU_STAMNES_H

#include "optics_utils.h"


/* @brief Hu and Stamnes cloud liquid water parameterization.*/
typedef struct HuStamnes
{
    double * a1; /*Exctinction coefficient parameter (radius, band).*/
    double * a2; /*Single-scatter albedo parameter (radius, band).*/
    double * a3; /*Asymmetry factor parameter (radius, band).*/
    double * band_limits; /*Lower/upper bounds of parameterization [cm-1] (band, 2).*/
    double * bands; /* Parameterization band centers [cm-1] (band).*/
    double * b1; /*Exctinction coefficient parameter (radius, band).*/
    double * b2; /*Single-scatter albedo parameter (radius, band).*/
    double * b3; /*Asymmetry factor parameter (radius, band).*/
    double * c1; /*Exctinction coefficient parameter (radius, band).*/
    double * c2; /*Single-scatter albedo parameter (radius, band).*/
    double * c3; /*Asymmetry factor parameter (radius, band).*/
    int last_ir_band; /*Index of last infrared (longwave) band.*/
    double max_radius; /*Maximum radius defined in parameterization [micron].*/
    double min_radius; /*Minimum radius defined in parameterization [micron].*/
    int num_bands; /*Number of bands.*/
    int num_radius_bins; /*Number of radius bins.*/
    double * radii; /*Radius bins [micron] for parameterization (radius).*/
} HuStamnes_t;


/* @brief Constructs a HuStamnes object.*/
void construct_liquid_optics(HuStamnes_t * self, char const * path);


/* @brief Destructs a HuStamnes object.*/
void destruct_liquid_optics(HuStamnes_t * self);


/* @brief Calculates cloud optics.*/
void calculate_liquid_optics(HuStamnes_t const self, double const water_concentration,
                             double const equivalent_radius,
                             OpticalProperties_t * optical_properties);


#endif
