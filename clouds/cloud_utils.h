#ifndef CLOUD_UTILS_H
#define CLOUD_UTILS_H


#include "debug.h"


/** @brief Create a map from wavenumber to cloud band.*/
int grid_band_mapping(
    int const grid_size, /**< Wavenumber grid size.*/
    fp_t const w0, /**< Wavenumber grid starting point [cm-1].*/
    fp_t const dw, /**< Wavenumber grid spacing [cm-1].*/
    int const num_bands, /**< Number of cloud bands.*/
    fp_t * const band_limits, /**< Cloud band limits [cm-1] (band, 2).*/
    int * mapping /**< Band number map (wavenumber).*/
);


#ifdef __NVCC__
__global__ void grid_band_mapping_d(
    int const grid_size, /**< Wavenumber grid size.*/
    fp_t const w0, /**< Wavenumber grid starting point [cm-1].*/
    fp_t const dw, /**< Wavenumber grid spacing [cm-1].*/
    int const num_bands, /**< Number of cloud bands.*/
    fp_t * const band_limits, /**< Cloud band limits [cm-1] (band, 2).*/
    int * mapping /**< Band number map (wavenumber).*/
);
#endif


/** @brief Expand optics from bands to the wavenumber grid.*/
int process_optics(
    int const grid_size, /**< Wavenumber grid size.*/
    int const num_layers, /**< Number of layers.*/
    int const num_bands, /**< Number of cloud bands.*/
    int const * mapping, /**< Band number map (wavenumber).*/
    fp_t const * tau, /**< Optical depth (band, layer).*/
    fp_t const * omega, /**< Single-scatter albedo (band, layer).*/
    fp_t const * g, /**< Asymmetry factor (band, layer).*/
    fp_t * optics_tau, /**< Optical depth (layer, wavenumber).*/
    fp_t * optics_omega, /**< Single-scatter albedo (layer, wavenumber).*/
    fp_t * optics_g /**< Asymmetry factor (layer, wavenumber).*/
);


#ifdef __NVCC__
__global__ void process_optics_d(
    int const grid_size, /**< Wavenumber grid size.*/
    int const num_layers, /**< Number of layers.*/
    int const num_bands, /**< Number of cloud bands.*/
    int const * mapping, /**< Band number map (wavenumber).*/
    fp_t const * tau, /**< Optical depth (band, layer).*/
    fp_t const * omega, /**< Single-scatter albedo (band, layer).*/
    fp_t const * g, /**< Asymmetry factor (band, layer).*/
    fp_t * optics_tau, /**< Optical depth (layer, wavenumber).*/
    fp_t * optics_omega, /**< Single-scatter albedo (layer, wavenumber).*/
    fp_t * optics_g /**< Asymmetry factor (layer, wavenumber).*/
);
#endif


#endif
