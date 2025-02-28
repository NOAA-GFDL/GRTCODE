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

#ifndef KERNELS_H_
#define KERNELS_H_

#include <stdint.h>
#include "grtcode_utilities.h"
#include "spectral_bin.h"


/** @brief Calculate pressure-shifted line center positions.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_line_centers(uint64_t const num_lines, /**< Number of molecular lines.*/
                      int const num_layers, /**< Number of atmospheric layers.*/
                      fp_t const * const v0, /**< Unshifted line center positions [cm-1] (line).*/
                      fp_t const * const delta, /**< Air-broadened pressure shift [cm-1 atm-1] (line).*/
                      fp_t const * const p, /**< Pressure [atm] (layer).*/
                      fp_t * const vnn /**< Pressure-shifted line center positions [cm-1] (layer, line).*/
                     );


/** @brief Calculate the total partition functions for each isotopologue.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_partition_functions(int const num_layers, /**< Number of atmospheric layers.*/
                             int const mol_id, /**< Molecule id.*/
                             int const num_iso, /**< Number of molecular isotopologues.*/
                             fp_t const * const t, /**< Temperature [K] (layer).*/
                             fp_t * const q /**< Total partition function (layer, isotopologue).*/
                            );


/** @brief Calculate temperature-corrected line intensities.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_line_strengths(uint64_t const num_lines, /**< Number of molecular lines.*/
                        int const num_layers, /**< Number of atmospheric layers.*/
                        int const num_iso, /**< Number of molecular isotopologues.*/
                        int const * const iso, /**< Isotopologue id (line).*/
                        fp_t const * const s0, /**< Uncorrected line strengths [cm] (line).*/
                        fp_t const * const vnn, /**< Line center position [cm-1] (line).*/
                        fp_t const * const en, /**< Lower state energies [cm-1] (line).*/
                        fp_t const * const t, /**< Temperature [K] (layer).*/
                        fp_t const * const q, /**< Total partition function (layer, isotopologue).*/
                        fp_t * const snn /**< Temperature-corrected line strengths [cm] (layer, line).*/
                       );


/** @brief Calculate lorentz halfwidths.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_lorentz_hw(uint64_t const num_lines, /**< Number of molecular lines.*/
                    int const num_layers, /**< Number of atmospheric layers.*/
                    fp_t const * const n, /**< Coefficient of temperature dependence of air-broadened halfwidths (line).*/
                    fp_t const * const yair, /**< Air-broadended halfwidths [cm-1 atm-1] at 296 K and 1 atm (line).*/
                    fp_t const * const yself, /**< Self-broadended halfwidths [cm-1 atm-1] at 296 K and 1 atm (line).*/
                    fp_t const * const t, /**< Temperature [K] (layer).*/
                    fp_t const * const p, /**< Pressure [atm] (layer).*/
                    fp_t const * const ps, /**< Partial pressure [atm] (layer).*/
                    fp_t * const gamma /**< Temperature and pressure corrected lorentz halfwidths [cm-1] (layer, line).*/
                   );


/** @brief Calculate doppler halfwidths.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_doppler_hw(uint64_t const num_lines, /**< Number of molecular lines.*/
                    int const num_layers, /**< Number of atmospheric layers.*/
                    fp_t const m, /**< Molecular mass [g].*/
                    fp_t const * const vnn, /**< Pressure-shifted line center position [cm-1] (layer, line).*/
                    fp_t const * const t, /**< Temperature [K] (layer).*/
                    fp_t * const alpha /**< Doppler halfwidths [cm-1] (layer, line).*/
                   );


/** @brief Sort the line parameters in order of line center wavenumber.
    @return GRTCODE_SUCCESS or an error code.*/
int sort_lines(uint64_t const num_lines, /**< Number of molecular lines.*/
               int const num_layers, /**< Number of molecular_lines.*/
               fp_t * const vnn, /**< Pressure-shifted line center positions [cm-1] (layer, line).*/
               fp_t * const snn, /**< Line strength [cm-1] (layer, line).*/
               fp_t * const gamma, /**< Lorentz halfwidth [cm-1] (layer, line).*/
               fp_t * const alpha /**< Doppler halfwidth [cm-1] (layer, line).*/
              );


/** @brief Calculate optical depths.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_optical_depth_bin_sweep(uint64_t const num_lines, /**< Number of molecular lines.*/
                                 int const num_layers, /**< Number of atmospheric layers.*/
                                 fp_t * const vnn, /**< Pressure-shifted line center positions [cm-1] (layer, line).*/
                                 fp_t * const snn, /**< Line strength [cm] (layer, line).*/
                                 fp_t * const gamma, /**< Lorentz halfwidth [cm-1] (layer, line).*/
                                 fp_t * const alpha, /**< Doppler halfwidth [cm-1] (layer, line).*/
                                 fp_t const * const n, /**< Integrated number density [cm-2] (layer).*/
                                 SpectralBins_t bins, /**< Spectral bins.*/
                                 fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                                );


/** @brief Calculate optical depths.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_optical_depth_line_sweep(uint64_t const num_lines, /**< Number of molecular lines.*/
                                  int const num_layers, /**< Number of atmospheric layers.*/
                                  fp_t * const vnn, /**< Pressure-shifted line center positions [cm-1] (layer, line).*/
                                  fp_t * const snn, /**< Line strength [cm] (layer, line).*/
                                  fp_t * const gamma, /**< Lorentz halfwidth [cm-1] (layer, line).*/
                                  fp_t * const alpha, /**< Doppler halfwidth [cm-1] (layer, line).*/
                                  fp_t const * const n, /**< Integrated number density [cm-2] (layer).*/
                                  SpectralBins_t bins, /**< Spectral bins.*/
                                  fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                                 );


/** @brief Calculate optical depths by sampling all lines.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_optical_depth_line_sample(uint64_t const num_lines, /**< Number of molecular lines.*/
                                   int const num_layers, /**< Number of atmospheric layers.*/
                                   fp_t * const vnn, /**< Pressure-shifted line center positions [cm-1] (layer, line).*/
                                   fp_t * const snn, /**< Line strength [cm] (layer, line).*/
                                   fp_t * const gamma, /**< Lorentz halfwidth [cm-1] (layer, line).*/
                                   fp_t * const alpha, /**< Doppler halfwidth [cm-1] (layer, line).*/
                                   fp_t const * const n, /**< Integrated number density [cm-2] (layer).*/
                                   SpectralBins_t const bins, /**< Spectral bins.*/
                                   fp_t * const tau, /**< Optical depth (layer, wavenumber).*/
                                   fp_t const * pedestal_lower_bound, /**< Pedestal lower bound [cm-1].*/
                                   fp_t const * pedestal_upper_bound /**< Pedestal upper bound [cm-1].*/
                                  );


/** @brief Calculate the optical depth contribution of the water vapor continuum.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_water_vapor_ctm_optical_depth(uint64_t const num_wpoints, /**< Spectral grid size.*/
                                       int const num_layers, /**< Number of atmospheric layers.*/
                                       fp_t * const tau, /**< Optical depth (layer, wavenumber).*/
                                       fp_t const * const CS, /**< (wavenumber).*/
                                       fp_t const * const T, /**< Temperature [K] (layer).*/
                                       fp_t const * const Ps, /**< Partial pressure [atm] (layer).*/
                                       fp_t const * const N, /**< Integrated number density [cm-2] (layer).*/
                                       fp_t const * const T0, /**< (wavenumber).*/
                                       fp_t const * const CF, /**< (wavenumber).*/
                                       fp_t const * const P, /**< Pressure [atm] (layer).*/
                                       fp_t const * const T0F /**< (wavenumber).*/
                                      );


/** @brief Calculate the optical depth contribution of the ozone continuum.
    @return GRTCODE_SUCCESS or an error code.*/
int calc_ozone_ctm_optical_depth(uint64_t const num_wpoints, /**< Spectral grid size.*/
                                 int const num_layers, /**< Number of atmospheric layers.*/
                                 fp_t const * const cross_section, /**< Cross section [cm2] (wavenumber).*/
                                 fp_t const * const N, /**< Integrated number density [cm-2] (layer).*/
                                 fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                                );


/** @brief Do a quadratic interpolation of line wing values in each bin
           (except for the last one.).
    @return GRTCODE_SUCCESS or an error code.*/
int interpolate(SpectralBins_t const bins, /**< Spectral bins.*/
                fp_t * const tau /**< Optical depth (layer, wavenumber).*/
               );


/** @brief Do a quadratic interpolation of line wing values in the last bin.
    @return GRTCODE_SUCCESS or an error code.*/
int interpolate_last_bin(SpectralBins_t const bins, /**< Spectral bins.*/
                         fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                        );


/** @brief Calculate the optical depth contribution of a CFC.
    @return SUCCESS or an error code.*/
int calc_cfc_optical_depth(uint64_t const num_wpoints, /**< Spectral grid size.*/
                           int const num_layers, /**< Number of atmospheric layers.*/
                           fp_t const * const n, /**< Integrated number density [cm^-2] (layers).*/
                           fp_t const * const x, /**< CFC abundance [ppmv] (levels).*/
                           fp_t const * const cross_section, /**< CFC cross section [cm^2] (wavenumber).*/
                           fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                          );


/** @brief Calculate the optical depth contribution of collision-induced absorption.
    @return SUCCESS or an error code.*/
int calc_cia_optical_depth(uint64_t const num_wpoints, /**< Spectral grid size.*/
                           int const num_layers, /**< Number of atmospheric layers.*/
                           fp_t const * const p, /**< Pressure [atm] (levels).*/
                           fp_t const * const t, /**< Temperature [K] (layers).*/
                           fp_t const * const x1, /**< Abundance [ppmv] of species one (levels).*/
                           fp_t const * const x2, /**< Abundance [ppmv] of species two (levels).*/
                           fp_t const * const cross_section, /**< Collision-induced absorption cross section [cm^4] (wavenumber).*/
                           fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                          );


#endif
