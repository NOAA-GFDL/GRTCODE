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

#include <string.h>
#include "cfcs.h"
#include "collision_induced_absorption.h"
#ifdef __NVCC__
#include "cuda_kernels.cuh"
#endif
#include "curtis_godson.h"
#include "debug.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "kernels.h"
#include "launch.h"
#include "molecules.h"
#include "molecules-internal.h"
#include "ozone_continuum.h"
#include "spectral_bin.h"
#include "water_vapor_continuum.h"
#include "water_vapor_continuum-internal.h"


/*Driver for optical depth calculation.*/
int launch(GasOptics_t * const gas_optics, fp_t *p, fp_t *t, fp_t * const tau)
{
    not_null(gas_optics);
    not_null(p);
    not_null(t);
    not_null(tau);

    /*Set pointers to input data.*/
    if (gas_optics->device == HOST_ONLY)
    {
        gas_optics->p = p;
        gas_optics->t = t;
        gas_optics->tau = tau;
    }
    else
    {
        gmemcpy(gas_optics->p, p, gas_optics->num_levels, gas_optics->device, FROM_HOST);
        gmemcpy(gas_optics->t, t, gas_optics->num_levels, gas_optics->device, FROM_HOST);
    }

    /*Zero out buffers used to accumulate results.*/
    gmemset(gas_optics->tau, 0 ,gas_optics->num_layers*gas_optics->grid.n,
            gas_optics->device);
    gmemset(gas_optics->bins.tau, 0, gas_optics->bins.isize*gas_optics->bins.num_layers,
            gas_optics->bins.device);

    /*Calculate the total number density of air molecules integrated across
      each layer.*/
    glaunch(calc_number_densities, gas_optics->num_layers, gas_optics->device,
            gas_optics->num_layers, gas_optics->p, gas_optics->n);

    /*Calculate integrated average layer quantities.*/
    glaunch(calc_pressures_and_temperatures, gas_optics->num_layers, gas_optics->device,
            gas_optics->num_layers, gas_optics->p, gas_optics->t, gas_optics->pavg,
            gas_optics->tavg);

    /*Loop over the molecules and calculate the optical depths.*/
    int m;
    for (m=0; m<gas_optics->num_molecules; ++m)
    {
        Molecule_t *mol = &(gas_optics->mols[m]);
        int index;
        catch(molecule_hash(mol->id, &index));
        char const *mesg = "Calculating spectra for %s.";
        log_info(mesg, mol->name);

        fp_t const plower = 0.;
        fp_t const pupper = 2.e4;
        fp_t const * pedestal_lower_bound = NULL;
        fp_t const * pedestal_upper_bound = NULL;
        if (mol->id == H2O)
        {
            pedestal_lower_bound = &plower;
            pedestal_upper_bound = &pupper;
        }

        /*Calculate the integrated average layer partial pressure.*/
        fp_t const *xp = &(gas_optics->x[index*gas_optics->num_levels]);
        glaunch(calc_partial_pressures_and_number_densities, gas_optics->num_layers,
                gas_optics->device, gas_optics->num_layers, gas_optics->p, xp,
                gas_optics->n, gas_optics->psavg, gas_optics->ns);

        /*Calculate pressure shifted line center positions.*/
        glaunch(calc_line_centers, mol->line_params.num_lines, gas_optics->device,
                mol->line_params.num_lines, gas_optics->num_layers, mol->line_params.vnn,
                mol->line_params.d, gas_optics->pavg, gas_optics->linecenter);

        /*Calculate total partition functions.*/
        glaunch(calc_partition_functions, mol->num_isotopologues, gas_optics->device,
                gas_optics->num_layers, mol->id, mol->num_isotopologues, gas_optics->tavg,
                mol->q);

        /*Calculate temperature-corrected line strengths.*/
        glaunch(calc_line_strengths, mol->line_params.num_lines, gas_optics->device,
                mol->line_params.num_lines, gas_optics->num_layers, mol->num_isotopologues,
                mol->line_params.iso, mol->line_params.snn, mol->line_params.vnn,
                mol->line_params.en, gas_optics->tavg, mol->q, gas_optics->snn);

        /*Calcluate temperature and pressure corrected lorentz half-widths.*/
        glaunch(calc_lorentz_hw, mol->line_params.num_lines, gas_optics->device,
                mol->line_params.num_lines, gas_optics->num_layers, mol->line_params.n,
                mol->line_params.yair, mol->line_params.yself, gas_optics->tavg,
                gas_optics->pavg, gas_optics->psavg, gas_optics->gamma);

        /*Calculate doppler half-widths.*/
        glaunch(calc_doppler_hw, mol->line_params.num_lines, gas_optics->device,
                mol->line_params.num_lines, gas_optics->num_layers, mol->mass,
                gas_optics->linecenter, gas_optics->tavg, gas_optics->alpha);

        /*Calculate the molecule's optical depths and add them to existing
          values.*/
        switch (gas_optics->optical_depth_method)
        {
            case wavenumber_sweep:
                glaunch(sort_lines, gas_optics->num_layers, gas_optics->device,
                        mol->line_params.num_lines, gas_optics->num_layers,
                        gas_optics->linecenter, gas_optics->snn, gas_optics->gamma,
                        gas_optics->alpha);
                glaunch(calc_optical_depth_bin_sweep, gas_optics->bins.n, gas_optics->device,
                        mol->line_params.num_lines, gas_optics->num_layers,
                        gas_optics->linecenter, gas_optics->snn, gas_optics->gamma,
                        gas_optics->alpha, gas_optics->ns, gas_optics->bins, gas_optics->tau);
                break;
            case line_sweep:
                glaunch(calc_optical_depth_line_sweep, mol->line_params.num_lines,
                        gas_optics->device, mol->line_params.num_lines, gas_optics->num_layers,
                        gas_optics->linecenter, gas_optics->snn, gas_optics->gamma,
                        gas_optics->alpha, gas_optics->ns, gas_optics->bins, gas_optics->tau);
                break;
            case line_sample:
                glaunch(calc_optical_depth_line_sample, mol->line_params.num_lines,
                        gas_optics->device, mol->line_params.num_lines, gas_optics->num_layers,
                        gas_optics->linecenter, gas_optics->snn, gas_optics->gamma,
                        gas_optics->alpha, gas_optics->ns, gas_optics->bins,
                        gas_optics->tau, pedestal_lower_bound, pedestal_upper_bound);
                break;
        }

        if (gas_optics->use_h2o_ctm && mol->id == H2O)
        {
            /*Calculate the water vapor continuum optical depths.*/
            glaunch(calc_water_vapor_ctm_optical_depth, gas_optics->bins.num_wpoints,
                    gas_optics->device, gas_optics->bins.num_wpoints, gas_optics->num_layers,
                    gas_optics->tau, gas_optics->h2o_cc.coefs[MTCKD25_S296], gas_optics->tavg,
                    gas_optics->psavg, gas_optics->ns, gas_optics->h2o_cc.coefs[CKDS],
                    gas_optics->h2o_cc.coefs[MTCKD25_F296], gas_optics->pavg,
                    gas_optics->h2o_cc.coefs[CKDF]);
        }
        else if (gas_optics->use_o3_ctm && mol->id == O3)
        {
            /*Calculate the ozone continuum optical depths.*/
            glaunch(calc_ozone_ctm_optical_depth, gas_optics->bins.num_wpoints,
                    gas_optics->device, gas_optics->bins.num_wpoints, gas_optics->num_layers,
                    gas_optics->o3_cc.cross_section, gas_optics->ns, gas_optics->tau);
        }
    }

    for (m=0; m<gas_optics->num_cfcs; ++m)
    {
        CfcCrossSection_t *cfc = &(gas_optics->cfcs[m]);
        int index = cfc->id;
        fp_t const *xp = &(gas_optics->x_cfc[index*gas_optics->num_levels]);
        char const *mesg = "Calculating spectra for %s.";
        log_info(mesg, cfc->name);

        /*Calculate CFC optical depths.*/
        glaunch(calc_cfc_optical_depth, gas_optics->bins.num_wpoints, gas_optics->device,
                gas_optics->bins.num_wpoints, gas_optics->num_layers, gas_optics->n, xp,
                cfc->cross_section, gas_optics->tau);
    }

    for (m=0; m<gas_optics->num_cias; ++m)
    {
        CollisionInducedAbsorption_t *cia = &(gas_optics->cia[m]);
        int index1 = cia->id[0];
        fp_t const *xp1 = &(gas_optics->x_cia[index1*gas_optics->num_levels]);
        int index2 = cia->id[1];
        fp_t const *xp2 = &(gas_optics->x_cia[index2*gas_optics->num_levels]);
        char const *mesg = "Calculating CIA spectra for %s - %s.";
        log_info(mesg, cia->name[0], cia->name[1]);

        /*Calculate collision-induced absorption optical depths.*/
        glaunch(calc_cia_optical_depth, gas_optics->bins.num_wpoints, gas_optics->device,
                gas_optics->bins.num_wpoints, gas_optics->num_layers, gas_optics->p,
                gas_optics->tavg, xp1, xp2, cia->cross_section, gas_optics->tau);
    }

    if (gas_optics->optical_depth_method != line_sample)
    {
        /*Interpolate line wing optical depth contributions.*/
        glaunch(interpolate, gas_optics->bins.n-1, gas_optics->device, gas_optics->bins,
                gas_optics->tau);
        glaunch(interpolate_last_bin, gas_optics->num_layers, gas_optics->device,
                gas_optics->bins, gas_optics->tau);
    }

    if (gas_optics->device != HOST_ONLY)
    {
        gmemcpy(tau, gas_optics->tau, gas_optics->num_layers*gas_optics->grid.n,
                gas_optics->device, FROM_DEVICE);
    }
    return GRTCODE_SUCCESS;
}
