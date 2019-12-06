/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef ATMOSPHERE_H_
#define ATMOSPHERE_H_

#include "grtcode_utilities.h"


/** @brief Atmospheric properties.*/
typedef struct Atmosphere
{
    int clean; /**< Assume the sky is clean.*/
    int clear; /**< Assume the sky is clear.*/
    fp_t alpha; /**< Surface albedo.*/
    int z; /**< Level lower bound index.*/
    int Z; /**< Level upper bound index.*/
    int num_levels; /**< Number of atmopsheric levels.*/
    int num_layers; /**< Number of atmophseric layers.*/
    SpectralGrid_t grid; /**< Spectral grid.*/
    int num_molecules; /**< Number of molecules.*/
    int num_cfcs; /**< Number of CFCs.*/
    int num_cias; /**< Number of CIAs.*/
    fp_t *level_pressure; /**< Pressure [atm] (level).*/
    fp_t *layer_pressure; /**< Pressure [atm] (layer).*/
    fp_t *level_temperature; /**< Temperature [K] (level).*/
    fp_t *layer_temperature; /**< Temperature [K] (layer).*/
    fp_t surface_temperature; /**< Surface temperature [K].*/
    fp_t total_solar_irradiance; /**< Total solar irradiance at TOA [W/m^2].*/
    fp_t solar_zenith_angle; /**< Cosine of solar zenith angle.*/
    fp_t *surface_albedo; /**< Surface albedo (wavenumber).*/
    fp_t *surface_emissivity; /**< Surface emissivity (wavenumber).*/
    fp_t **ppmv; /**< Molecular abundance [ppmv] (molecule, level).*/
    fp_t **cfc_ppmv; /**< CFC abundance [ppmv] (CFC, level).*/
    fp_t **cia_ppmv; /**< CIA abindance [ppmv] (molecule, level).*/
    fp_t *aerosol_optical_depth; /**< Aerosol optical depth (layer, wavenumber).*/
    fp_t *aerosol_single_scatter_albedo; /**< Aerosol single-scatter albedo (layer, wavenumber).*/
    fp_t *aerosol_asymmetry_factor; /**< Aerosol asymmetry factory (layer, wavenumber).*/
    fp_t *liquid_water_path; /**< Liquid water path [g/m^2] (layer).*/
    fp_t *liquid_water_droplet_radius; /**< Liquid water equivalent radiu [microns] (layer).*/
} Atmosphere_t;


/**@ brief Dimension indices.*/
enum dimensions
{
    LEVEL = 0,
    NUM_DIMS
};


/** @brief Flux variable indices.*/
enum variables
{
    RLU = 0,
    RLD,
    RSU,
    RSD,
    NUM_VARS
};


/** @brief Output file object.*/
typedef struct Output
{
    int ncid;
    int dimid[NUM_DIMS];
    int varid[NUM_VARS];
} Output_t;


/** @brief Container for CFC/HFC arguments.*/
typedef struct Cfc
{
    int id;
    char *path;
} Cfc_t;


/**@ brief Reserve memory and read in atmospheric data.*/
void create_atmosphere(Atmosphere_t * const atm, /**< Atmosphere object.*/
                       char const * const filepath, /**< Input data file.*/
                       int const * const molecules, /**< Array of molecule ids.*/
                       int const num_molecules, /**< Number of molecules.*/
                       Cfc_t const * const cfc, /**< Array of Cfc_t objects.*/
                       int const num_cfcs, /**< Number of CFCs.*/
                       int const * const cias, /**< Array of CIA ids.*/
                       int const num_cias /**< NUmber of CIAs.*/
                      );


/** @brief Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm /*Atmosphere object.*/
                       );


/** @brief Create an output file and write metadata.*/
Output_t create_flux_file(char const * const filepath, /**< File path.*/
                          Atmosphere_t const * const atm /**< Atmosphere object.*/
                         );


/** @brief Close output file.*/
void close_flux_file(Output_t const * const o /**< Output object.*/
                    );


/** @brief Write a column of fluxes to the output file.*/
void write_fluxes(Output_t const * const o, /**< Output object.*/
                  int const index, /**< Variable index.*/
                  fp_t const * const flux /**< Fluxes [W/m^2] (level).*/
                 );


#endif
