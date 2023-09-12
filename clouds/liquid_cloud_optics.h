#ifndef HU_STAMNES_H
#define HU_STAMNES_H


#include "debug.h"
#include "device.h"


/* @brief Liquid water cloud optics parameterization.*/
typedef struct LiquidCloudOptics
{
    fp_t * a1; /*Exctinction coefficient parameter (radius, band).*/
    fp_t * a2; /*Single-scatter albedo parameter (radius, band).*/
    fp_t * a3; /*Asymmetry factor parameter (radius, band).*/
    fp_t * band_limits; /*Lower/upper bounds of parameterization [cm-1] (band, 2).*/
    fp_t * bands; /* Parameterization band centers [cm-1] (band).*/
    fp_t * b1; /*Exctinction coefficient parameter (radius, band).*/
    fp_t * b2; /*Single-scatter albedo parameter (radius, band).*/
    fp_t * b3; /*Asymmetry factor parameter (radius, band).*/
    fp_t * c1; /*Exctinction coefficient parameter (radius, band).*/
    fp_t * c2; /*Single-scatter albedo parameter (radius, band).*/
    fp_t * c3; /*Asymmetry factor parameter (radius, band).*/
    Device_t device; /*Device that the code is run on.*/
    int last_ir_band; /*Index of last infrared (longwave) band.*/
    fp_t max_radius; /*Maximum radius defined in parameterization [micron].*/
    fp_t min_radius; /*Minimum radius defined in parameterization [micron].*/
    int num_bands; /*Number of bands.*/
    int num_radius_bins; /*Number of radius bins.*/
    fp_t * radii; /*Radius bins [micron] for parameterization (radius).*/
} LiquidCloudOptics_t;


/* @brief Constructs a HuStamnes object.*/
EXTERN int create_liquid_optics(
    LiquidCloudOptics_t * self, /**< Liquid optics object.*/
    char const * path, /**< Path to liquid optics input file.*/
    Device_t const * const device /**< Device to run on.*/
);


/* @brief Destructs a HuStamnes object.*/
EXTERN int destroy_liquid_optics(
    LiquidCloudOptics_t * self /**< Liquid optics object.*/
);


/* @brief Calculates cloud optics.*/
int calculate_liquid_optics(
    fp_t const min_radius, /**< Minimum radius [micron].*/
    fp_t const max_radius, /**< Maximum radius [micron].*/
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_bands, /**< Number of bands.*/
    fp_t const * a1, /**< a1 parameter (radius, band).*/
    fp_t const * b1, /**< b1 parameter (radisu, band).*/
    fp_t const * c1, /**< c1 parameter (radius, band).*/
    fp_t const * a2, /**< a2 parameter (radius, band).*/
    fp_t const * b2, /**< b2 parameter (radius, band).*/
    fp_t const * c2, /**< c2 parameter (radius, band).*/
    fp_t const * a3, /**< a3 parameter (radius, band).*/
    fp_t const * b3, /**< b3 parameter (radius, band).*/
    fp_t const * c3, /**< c3 parameter (radius, band).*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * water_concentration, /**< Water concentration [g m-3] (band, layer).*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const * thickness, /**< Layer thickness [m] (layer).*/
    fp_t * optical_depth, /**< Optical depth (band, layer).*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo (band, layer).*/
    fp_t * asymmetry_factor /**< Asymmetry factor (band, layer).*/
);


#ifdef __NVCC__
/* @brief Calculates liquid cloud optics for all bands.*/
__global__ void calculate_liquid_optics_d(
    fp_t const min_radius, /**< Minimum radius [micron].*/
    fp_t const max_radius, /**< Maximum radius [micron].*/
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_bands, /**< Number of bands.*/
    fp_t const * a1, /**< a1 parameter (radius, band).*/
    fp_t const * b1, /**< b1 parameter (radisu, band).*/
    fp_t const * c1, /**< c1 parameter (radius, band).*/
    fp_t const * a2, /**< a2 parameter (radius, band).*/
    fp_t const * b2, /**< b2 parameter (radius, band).*/
    fp_t const * c2, /**< c2 parameter (radius, band).*/
    fp_t const * a3, /**< a3 parameter (radius, band).*/
    fp_t const * b3, /**< b3 parameter (radius, band).*/
    fp_t const * c3, /**< c3 parameter (radius, band).*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * water_concentration, /**< Water concentration [g m-3] (band, layer).*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const * thickness, /**< Layer thickness [m] (layer).*/
    fp_t * optical_depth, /**< Optical depth (band, layer).*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo (band, layer).*/
    fp_t * asymmetry_factor /**< Asymmetry factor (band, layer).*/
);
#endif


#endif
