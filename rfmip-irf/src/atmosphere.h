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
    int x; /**< Column lower bound index.*/
    int z; /**< Level lower bound index.*/
    int num_columns; /**< Number of atmospheric columns.*/
    int num_levels; /**< Number of atmopsheric levels.*/
    int num_layers; /**< Number of atmophseric layers.*/
    uint64_t num_lw_wavenumber; /**< Number of wavenumber points for longwave solver.*/
    uint64_t num_sw_wavenumber; /**< Number of wavenumber points for shortwave solver.*/
    int num_molecules; /**< Number of molecules.*/
    int num_cfcs; /**< Number of CFCs.*/
    int num_cias; /**< Number of CIAs.*/
    fp_t *level_pressure; /**< Pressure [atm] (column, level).*/
    fp_t *layer_pressure; /**< Pressure [atm] (column, layer).*/
    fp_t *level_temperature; /**< Temperature [K] (column, level).*/
    fp_t *layer_temperature; /**< Temperature [K] (column, layer).*/
    fp_t *surface_temperature; /**< Surface temperature [K] (column).*/
    fp_t *total_solar_irradiance; /**< Total solar irradiance at TOA [W/m^2] (column).*/
    fp_t *solar_zenith_angle; /**< Cosine of solar zenith angle (column).*/
    fp_t *surface_albedo; /**< Surface albedo (column, wavenumber).*/
    fp_t *surface_emissivity; /**< Surface emissivity (column, wavenumber).*/
    fp_t **ppmv; /**< Molecular abundance [ppmv] (molecule, column, level).*/
    fp_t **cfc_ppmv; /**< CFC abundance [ppmv] (CFC, column, level).*/
    fp_t **cia_ppmv; /**< CIA abindance [ppmv] (molecule, column, level).*/
} Atmosphere_t;


/**@ brief Dimension indices.*/
enum dimensions
{
    EXPT = 0,
    SITE,
    LEVEL,
    LAYER,
    TIME,
    NUM_DIMS
};


/** @brief Flux variable indices.*/
enum variables
{
    RLU = 0,
    RLD,
    RSU,
    RSD,
    RLHR,
    RSHR,
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
    int use_equivalent_ppmv;
} Cfc_t;


/**@ brief Reserve memory and read in atmospheric data.*/
void create_atmosphere(Atmosphere_t * const atm, /**< Atmosphere object.*/
                       char const * const filepath, /**< Input data file.*/
                       int const experiment, /**< Experiment index.*/
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
void write_output(Output_t const * const o, /**< Output object.*/
                  int const index, /**< Variable index.*/
                  fp_t const * const data, /**< Data.*/
                  size_t const * const start, /**< Corner.*/
                  size_t const * const count /**< Edge lengths.*/
                 );


void fixed_dynamic_heating(char const * const path, fp_t * const dq);


#endif
