/** @file*/
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

/** @file */
#ifndef GAS_OPTICS_H_
#define GAS_OPTICS_H_

#include <stdint.h>
#include "cfcs.h"
#include "collision_induced_absorption.h"
#include "grtcode_utilities.h"
#include "molecules.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "water_vapor_continuum.h"


/** @brief Maximum file path size.*/
#define DIR_PATH_LEN 1024


/** @brief Gas optics.*/
typedef struct GasOptics
{
    Device_t device; /**< Id of the device associated with this object.*/
    int num_levels; /**< Number of atmospheric levels.*/
    int num_layers; /**< Number of atmospheric layers (= num_levels - 1).*/

    int num_molecules; /**< Number of molecules.*/
    uint64_t molecule_bit_field; /**< Bit field used to determine which molecules are currently in use.*/
    Molecule_t mols[NUM_MOLS]; /**< Array of molecule structures.*/

    int num_cfcs; /**< Number of cfcs.*/
    uint64_t cfc_bit_field; /**< Bit field used to determine which cfcs are currently in use.*/
    CfcCrossSection_t cfcs[NUM_CFCS]; /**< CFC cross section data structures.*/
    fp_t *x_cfc; /**< CFC abundance (CFC, level).*/

    int num_cias; /**< Number of collision-induced absorption continua.*/
    uint64_t cia_bit_field; /**< Bit field used to determine with species are currently in use.*/
    CollisionInducedAbsorption_t cia[MAX_NUM_CIAS]; /**< Collision-induce absorption structures.*/
    fp_t *x_cia; /**< Collision-induced absorption abundances (molecule, level).*/

    char h2o_ctm_dir[DIR_PATH_LEN]; /**< Path to the water vapor continuum directory.*/
    int use_h2o_ctm; /**< Flag indicating if using the water vapor continuum is used.*/
    WaterVaporContinuumCoefs_t h2o_cc; /**< Water vapor continuum coefficients.*/

    char o3_ctm_file[DIR_PATH_LEN]; /**< Path to the ozone continuum file.*/
    int use_o3_ctm; /**< Flag indicating if the ozone continuum is used.*/
    OzoneContinuumCoefs_t o3_cc; /**< Oone continuum coefficients.*/

    SpectralGrid_t grid; /**< Spectral grid.*/
    SpectralBins_t bins; /**< Spectral bins.*/
    char hitran_path[DIR_PATH_LEN]; /**< Path to the HITRAN database file.*/
    double wcutoff; /**< Cutoff from spectral line center [cm-1].*/
    int optical_depth_method; /**< Flag specifying which method will be used to calculate the optical depths.*/
    fp_t *x; /**< Abundance (molecule, level).*/
    fp_t *n; /**< Integrated number density [cm-2] (layer).*/
    fp_t *pavg; /**< Pressure [atm] (layer).*/
    fp_t *tavg; /**< Temperature [K] (layer).*/
    fp_t *psavg; /**< Molecular partial pressure [atm] (layer).*/
    fp_t *ns; /**< Integrated number density [cm-2] of a particular species (layer).*/
    fp_t *linecenter; /**< Pressure-shifted line center position [cm-1] (layer, line).*/
    fp_t *snn; /**< Line strength [cm-1] (layer, line).*/
    fp_t *gamma; /**< Temperature- and pressure-corrected lorentz half-width [cm-1] (layer, line).*/
    fp_t *alpha; /**< Doppler half-width [cm-1] (layer, line).*/
    fp_t *p; /**< Pressure [atm] (level).*/
    fp_t *t; /**< Temperature [K] (level).*/
    fp_t *tau; /**< Optical depth (layer, wavenumber).*/
} GasOptics_t;


/** @brief Flags used to specifiy which method is used to calculate the optical depths.*/
enum OpticalDepthMethod
{
    wavenumber_sweep, /**< RFM-like sweep through spectral bins.*/
    line_sweep, /**< RFM-like sweep, but through spectral lines.*/
    line_sample /**< Sampling of each spectral line at the desired resolution.*/
};


/** @brief Reserve memory for molecular lines.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int create_gas_optics(GasOptics_t * const gas_optics, /**< Gas optics.*/
                             int const num_levels, /**< Number of atmospheric levels.*/
                             SpectralGrid_t const * const grid, /**< Spectral grid.*/
                             Device_t const * const device, /**< Device.*/
                             char const * const hitran_path, /**< Path to HITRAN database file.*/
                             char const * const h2o_ctm_dir, /**< Path to water vapor continuum directory.*/
                             char const * const o3_ctm_file, /**< Path to ozone continuum file.*/
                             double const * const wcutoff, /**< Cutoff from line center [1/cm].*/
                             int const * const optical_depth_method /**< Method to use to calculate optical depths.*/
                            );


/** @brief Free memory for the molecular lines.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int destroy_gas_optics(GasOptics_t * const gas_optics /**< Gas optics.*/
                             );


/** @brief Add a molecule.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int add_molecule(GasOptics_t * const gas_optics, /**< Gas optics.*/
                        int const molecule_id, /**< Molecule id.*/
                        double const * const min_line_center, /**< Lower bound [cm-1] for spectral line centers.*/
                        double const * const max_line_center /**< Upper bound [cm-1] for spectral line centers.*/
                       );


/** @brief Update a molecule's ppmv.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int set_molecule_ppmv(GasOptics_t * const gas_optics, /**< Gas optics.*/
                             int const molecule_id, /**< Molecule id.*/
                             fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                            );


/** @brief Add a CFC.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int add_cfc(GasOptics_t * const gas_optics, /**< Gas optics.*/
                   int const cfc_id, /**< CFC id.*/
                   char const * const filepath /**< Path to CFC cross section csv file.*/
                  );


/** @brief Update a CFC's ppmv.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int set_cfc_ppmv(GasOptics_t * const gas_optics, /**< Gas optics.*/
                        int const cfc_id, /**< CFC id.*/
                        fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                       );


/** @brief Activate collision-induced absorption between two species.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int add_cia(GasOptics_t * const gas_optics, /**< Gas optics.*/
                   int const species1, /**< Id of species.*/
                   int const species2, /**< Id of species.*/
                   char const * const filepath /**< Path to cross section csv file.*/
                  );


/** @brief Update a CIA species' ppmv.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int set_cia_ppmv(GasOptics_t * const ml, /**< Gas optics.*/
                        int const cia_id, /**< CIA species id.*/
                        fp_t const * const ppmv /**< Abundance [ppmv] (level).*/
                       );


/** @brief Calcluate the total optical depth in each layer at each spectral grid point.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int calculate_optical_depth(GasOptics_t * const gas_optics, /**< Gas optics.*/
                                   fp_t * const pressure, /**< Pressure [mb] (level).*/
                                   fp_t * const temperature, /**< Temperature [K] (level).*/
                                   Optics_t * const optics /**< Optics object.*/
                                  );


/** @brief Get the number of molecules.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int get_num_molecules(GasOptics_t const * const gas_optics, /**< Gas optics.*/
                             int * const n /**< Number of molecules.*/
                            );


#endif
