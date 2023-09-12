#ifndef ICE_CLOUD_OPTICS_H
#define ICE_CLOUD_OPTICS_H


#include "debug.h"
#include "device.h"
#include "optics_utils.h"


/* @brief Ice water cloud optics parameterizations.*/
typedef struct IceCloudOptics
{
    fp_t * a; /*a parameters from equations 4a/5a (band, 6).*/
    fp_t * b; /*b parameters from equations 4b/5b (band, 6).*/
    fp_t * band_limits; /*Parameterization band limits [cm-1] (bands, 2).*/
    fp_t * bands; /*Parameterization band centers [cm-1] (band).*/
    fp_t * c; /*c parameters from equations 4c/5c (radius, band, 6).*/
    Device_t device; /*Device to run the code on.*/
    int last_ir_band; /*Index of last infrared (longwave) band.*/
    int num_bands; /*Number or bands.*/
    int num_order; /*Number of polynomial orders.*/
    int num_radius_bins; /*Number of radius bins.*/
    fp_t * radii; /*Radius bins [micron] for the parameterization (radius, 2).*/
} IceCloudOptics_t;


/* @brief Constructs a IceCloudOptics object.*/
EXTERN int create_ice_optics(
    IceCloudOptics_t * self, /**< Ice optics object.*/
    char const * path, /**< Path to the ice optics input file.*/
    Device_t const * const device /**< Device to run on.*/
);


/* @brief Constructs a IceCloudOptics object.*/
EXTERN int destroy_ice_optics(
    IceCloudOptics_t * self /**< Ice optics object.*/
);


/* @brief Calculates cloud optics.*/
int calculate_ice_optics(
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_order, /**< Number of polynomial orders.*/
    int const num_bands, /**< Number of bands.*/
    int const last_ir_band, /**< Index of the last infrared band.*/
    fp_t const * a, /**< a parameter (band, order).*/
    fp_t const * b, /**< b parameter (band, order).*/
    fp_t const * c, /**< c parameter (band, order).*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * ice_concentration, /**< Ice concentration [g m-3] (band, layer).*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const scale_factor, /**< Scaling factor.*/
    fp_t const * temperature, /**< Temperature [K] (layer).*/
    fp_t const * thickness, /**< Layer thickness [m] (layer).*/
    fp_t * optical_depth, /**< Optical depth (band, layer).*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo (band, layer).*/
    fp_t * asymmetry_factor /**< Asymmetry factor (band, layer).*/
);


#ifdef __NVCC__
/* @brief Calculates cloud optics.*/
__global__ void calculate_ice_optics_d(
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_order, /**< Number of polynomial orders.*/
    int const num_bands, /**< Number of bands.*/
    int const last_ir_band, /**< Index of the last infrared band.*/
    fp_t const * a, /**< a parameter (band, order).*/
    fp_t const * b, /**< b parameter (band, order).*/
    fp_t const * c, /**< c parameter (band, order).*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * ice_concentration, /**< Ice concentration [g m-3] (band, layer).*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const scale_factor, /**< Scaling factor.*/
    fp_t const * temperature, /**< Temperature [K] (layer).*/
    fp_t const * thickness, /**< Layer thickness [m] (layer).*/
    fp_t * optical_depth, /**< Optical depth (band, layer).*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo (band, layer).*/
    fp_t * asymmetry_factor /**< Asymmetry factor (band, layer).*/
);
#endif


#endif
