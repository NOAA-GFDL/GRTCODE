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
#ifdef __NVCC__

#ifndef CUDA_KERNELS_CUH_
#define CUDA_KERNELS_CUH_

#include <stdint.h>
#include "floating_point_type.h"
#include "spectral_bin.h"


/** @brief Calculate pressure-shifted line center positions.*/
__global__ void calc_line_centers_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                    int const num_layers, /*Number of atmospheric layers.*/
                                    fp_t const * const v0, /*Unshifted line center
                                                             positions [1/cm] (lines).*/
                                    fp_t const * const delta, /*Air-broadened pressure
                                                                shift [1/(cm*atm)] (lines).*/
                                    fp_t const * const p, /*Pressure [atm] (layers).*/
                                    fp_t * const vnn /*Pressure-shifted line center
                                                       positions [1/cm] (layers,lines).*/
                                   );


/** @brief Calculate the total partition functions for each isotopologue.*/
__global__ void calc_partition_functions_d(int const num_layers, /*Number of atmospheric layers.*/
                                           int const mol_id, /*Molecule id.*/
                                           int const num_iso, /*Number of molecular isotopologues.*/
                                           fp_t const * const t, /*Temperature [K] (layers).*/
                                           fp_t * const q /*Total partition function.*/
                                          );


/** @brief Calculate temperature-corrected line intensities.*/
__global__ void calc_line_strengths_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                      int const num_layers, /*Number of atmospheric layers.*/
                                      int const num_iso, /*Number of molecular
                                                           isotopologues.*/
                                      int const * const iso, /*Isotopologue id (lines).*/
                                      fp_t const * const s0, /*Uncorrected line strengths
                                                               [1/cm] (lines).*/
                                      fp_t const * const vnn, /*Line center position [1/cm]
                                                                (lines).*/
                                      fp_t const * const en, /*Lower state energies [1/cm]
                                                               (lines).*/
                                      fp_t const * const t, /*Temperature [K] (layers).*/
                                      fp_t const * const q, /*Total partition function.*/
                                      fp_t * const snn /*Temperature-corrected line
                                                         strengths [1/cm] (layers,lines).*/
                                     );


/** @brief Calculate lorentz halfwidths.*/
__global__ void calc_lorentz_hw_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                  int const num_layers, /*Number of atmospheric layers.*/
                                  fp_t const * const n, /*Coefficient of temperature
                                                          dependence of air-broadened
                                                          halfwidths (lines).*/
                                  fp_t const * const yair, /*Air-broadended halfwidths
                                                             [1/(cm*atm)] at 296K and 1atm.*/
                                  fp_t const * const yself, /*Self-broadended halfwidths
                                                             [1/(cm*atm)] at 296K and 1atm.*/
                                  fp_t const * const t, /*Temperature [K] (layers).*/
                                  fp_t const * const p, /*Pressure [atm] (layers).*/
                                  fp_t const * const ps, /*Partial pressure [atm] (layers).*/
                                  fp_t * const gamma /*Temperature and pressure corrected
                                                       lorentz halfwidths [1/cm]
                                                       (layers,lines).*/
                                 );


/** @brief Calculate doppler halfwidths.*/
__global__ void calc_doppler_hw_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                  int const num_layers, /*Number of atmospheric layers.*/
                                  fp_t const m, /*Molecular mass [g].*/
                                  fp_t const * const vnn, /*Pressure-shifted line center
                                                            position [1/cm] (layers,lines).*/
                                  fp_t const * const t, /*Temperature [K] (layers).*/
                                  fp_t * const alpha /*Doppler halfwidths [1/cm]
                                                       (layers,lines).*/
                                 );


/** @brief Sort the line parameters in order of line center wavenumber.*/
__global__ void sort_lines_d(uint64_t const num_lines, /*Number of molecular lines.*/
                             int const num_layers, /*Number of molecular_lines.*/
                             fp_t * const vnn, /*Pressure-shifted line
                                                 center positions [1/cm].
                                                 (layers,lines).*/
                             fp_t * const snn, /*Line strength [1/cm]
                                                 (layers,lines).*/
                             fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                   (layers,lines).*/
                             fp_t * const alpha /*Doppler halfwidth [1/cm]
                                                  (layers,lines).*/
                            );


/** @brief Calculate optical depths.*/
__global__ void calc_optical_depth_bin_sweep_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                               int const num_layers, /*Number of atmospheric layers.*/
                                               fp_t * const vnn, /*Pressure-shifted line
                                                                   center positions [1/cm].
                                                                   (layers,lines).*/
                                               fp_t * const snn, /*Line strength [1/cm]
                                                                   (layers,lines).*/
                                               fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                                     (layers,lines).*/
                                               fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                                     (layers,lines).*/
                                               fp_t const * const n, /*Integrated number density
                                                                       [cm^-2] (layers).*/
                                               SpectralBins_t bins, /*Spectral bins.*/
                                               fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                              );


/** @brief Calculate optical depths.*/
__global__ void calc_optical_depth_line_sweep_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                                int const num_layers, /*Number of atmospheric layers.*/
                                                fp_t * const vnn, /*Pressure-shifted line
                                                                    center positions [1/cm].
                                                                    (layers,lines).*/
                                                fp_t * const snn, /*Line strength [1/cm]
                                                                    (layers,lines).*/
                                                fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                                      (layers,lines).*/
                                                fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                                      (layers,lines).*/
                                                fp_t const * const n, /*Integrated number density
                                                                        [cm^-2] (layers).*/
                                                SpectralBins_t bins, /*Spectral bins.*/
                                                fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                               );


/** @brief Calculate optical depths by sampling all lines.*/
__global__ void calc_optical_depth_line_sample_d(uint64_t const num_lines, /*Number of molecular lines.*/
                                                 int const num_layers, /*Number of atmospheric layers.*/
                                                 fp_t * const vnn, /*Pressure-shifted line
                                                                     center positions [1/cm].
                                                                     (layers,lines).*/
                                                 fp_t * const snn, /*Line strength [1/cm]
                                                                     (layers,lines).*/
                                                 fp_t * const gamma, /*Lorentz halfwidth [1/cm]
                                                                       (layers,lines).*/
                                                 fp_t * const alpha, /*Doppler halfwidth [1/cm]
                                                                       (layers,lines).*/
                                                 fp_t const * const n, /*Integrated number density
                                                                         [cm^-2] (layers).*/
                                                 SpectralBins_t const bins, /*Spectral bins.*/
                                                 fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                                );


/** @brief Calculate the optical depth contribution of the water vapor
           continuum.*/
__global__ void calc_water_vapor_ctm_optical_depth_d(uint64_t const num_wpoints, /*Spectral grid size.*/
                                                     int const num_layers, /*Number of atmospheric layers.*/
                                                     fp_t * const tau, /*Optical depth (layer,wavenumber).*/
                                                     fp_t const * const CS,
                                                     fp_t const * const T,
                                                     fp_t const * const Ps,
                                                     fp_t const * const N,
                                                     fp_t const * const T0,
                                                     fp_t const * const CF,
                                                     fp_t const * const P,
                                                     fp_t const * const T0F
                                                    );


/** @brief Calculate the optical depth contribution of the ozone continuum.*/
__global__ void calc_ozone_ctm_optical_depth_d(uint64_t const num_wpoints, /*Spectral grid size.*/
                                               int const num_layers, /*Number of atmospheric layers.*/
                                               fp_t const * const cross_section,
                                               fp_t const * const N,
                                               fp_t * const tau /*Optical depth (layer,wavenumber).*/
                                              );


/** @brief Do a quadratic interpolation of line wing values in each bin
           (except for the last one.).*/
__global__ void interpolate_d(SpectralBins_t const bins,
                              fp_t * const tau
                             );


/** @brief Do a quadratic interpolation of line wing values in each bin
           (except for the last one.).*/
__global__ void interpolate_last_bin_d(SpectralBins_t const bins,
                                       fp_t * const tau
                                      );


/** @brief Calculate the optical depth contribution of a CFC.*/
__global__ void calc_cfc_optical_depth_d(uint64_t const num_wpoints, /**< Spectral grid size.*/
                                         int const num_layers, /**< Number of atmospheric layers.*/
                                         fp_t const * const n, /**< Integrated number density [cm^-2] (layers).*/
                                         fp_t const * const x, /**< CFC abundance [ppmv] (levels).*/
                                         fp_t const * const cross_section, /**< CFC cross section [cm^2] (wavenumber).*/
                                         fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                                        );


/** @brief Calculate the optical depth contribution of collision-induced absorption.*/
__global__ void calc_cia_optical_depth_d(uint64_t const num_wpoints, /**< Spectral grid size.*/
                                         int const num_layers, /**< Number of atmospheric layers.*/
                                         fp_t const * const p, /**< Pressure [atm] (levels).*/
                                         fp_t const * const t, /**< Temperature [K] (layers).*/
                                         fp_t const * const x1, /**< Abundance [ppmv] of species one (levels).*/
                                         fp_t const * const x2, /**< Abundance [ppmv] of species two (levels).*/
                                         fp_t const * const cross_section, /**< Collision-induced absorption cross section [cm^4] (wavenumber).*/
                                         fp_t * const tau /**< Optical depth (layer, wavenumber).*/
                                        );


#endif
#endif
