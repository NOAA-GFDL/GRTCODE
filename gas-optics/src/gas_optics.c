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

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cfcs.h"
#include "cfcs-internal.h"
#include "collision_induced_absorption.h"
#include "collision_induced_absorption-internal.h"
#include "debug.h"
#include "gas_optics.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "launch.h"
#include "molecules.h"
#include "molecules-internal.h"
#include "ozone_continuum.h"
#include "ozone_continuum-internal.h"
#include "spectral_bin.h"
#include "spectral_bin-internal.h"
#include "tips2017.h"
#include "water_vapor_continuum.h"
#include "water_vapor_continuum-internal.h"


static double const MIN_CUTOFF = 1.; /*Smallest cut-off [cm-1] from a line center allowed.*/
static double const MAX_CUTOFF = 50.; /*Larget cut-off [cm-1] from a line center allowed.*/
static int const MAX_NUM_LINES = 600000; /*Largest number of spectral lines per molecule allowed.*/
static double const DEFAULT_CUTOFF = 25.; /*Default cut-off [cm-1] from a line center.*/


/*Reserve memory for gas optics.*/
EXTERN int create_gas_optics(GasOptics_t * const gas_optics, int const num_levels,
                             SpectralGrid_t const * const grid,
                             Device_t const * const device,
                             char const * const hitran_path,
                             char const * const h2o_ctm_dir,
                             char const * const o3_ctm_file,
                             double const * const wcutoff,
                             int const * const optical_depth_method)
{
    not_null(gas_optics);
    in_range(num_levels, MIN_NUM_LEVELS, MAX_NUM_LEVELS);
    gas_optics->num_levels = num_levels;
    gas_optics->num_layers = num_levels - 1;
    char const *mesg = "Molecular lines properties:\n\tnumber of levels: %d\n\t"
                       "number of layers: %d";
    log_info(mesg, gas_optics->num_levels, gas_optics->num_layers);
    not_null(grid);
    gas_optics->grid = *grid;
    not_null(device);
    gas_optics->device = *device;

    /*Create the spectral bins.*/
    double bin_width = 1.;
    catch(create_spectral_bins(&(gas_optics->bins), gas_optics->num_layers, gas_optics->grid.w0,
                               gas_optics->grid.n, gas_optics->grid.dw, bin_width,
                               gas_optics->device));
    mesg = "Spectral bin properties:\n\tnumber of bins: %zu\n\t"
           "bin width: %e\n\tspectral grid points per bin: %d\n\t"
           "interpolation: %d\n\tspectral gid points in last bin:"
           " %d\n\tinterpolation in last bin: %d";
    log_info(mesg, gas_optics->bins.n, bin_width, gas_optics->bins.ppb,
             gas_optics->bins.do_interp, gas_optics->bins.last_ppb,
             gas_optics->bins.do_last_interp);

    /*Store the path to the hitran database file.*/
    not_null(hitran_path);
    snprintf(gas_optics->hitran_path, DIR_PATH_LEN, "%s", hitran_path);
    mesg = "Using HITRAN database file %s.";
    log_info(mesg, gas_optics->hitran_path);

    /*Set the molecular line cutoff.*/
    if (wcutoff != NULL)
    {
        in_range(*wcutoff, MIN_CUTOFF, MAX_CUTOFF);
        gas_optics->wcutoff = *wcutoff;
    }
    else
    {
        gas_optics->wcutoff = DEFAULT_CUTOFF;
    }
    mesg = "Using spectral line cut-off of %e [1/cm].";
    log_info(mesg, gas_optics->wcutoff);

    /*Set the method that will be used to calculate the optical depths.*/
    if (optical_depth_method != NULL)
    {
        in_range(*optical_depth_method, wavenumber_sweep, line_sample);
        gas_optics->optical_depth_method = *optical_depth_method;
    }
    else
    {
        gas_optics->optical_depth_method = wavenumber_sweep;
    }

    /*Prepare to add molecules/cfcs/cias.*/
    gas_optics->num_molecules = 0;
    gas_optics->molecule_bit_field = 0;
    gas_optics->num_cfcs = 0;
    gas_optics->cfc_bit_field = 0;
    gas_optics->num_cias = 0;
    gas_optics->cia_bit_field = 0;

    /*Pepare water vapor continuum.*/
    gas_optics->use_h2o_ctm = 0;
    if (h2o_ctm_dir != NULL)
    {
        if (strcmp(h2o_ctm_dir, "none") != 0)
        {
            gas_optics->use_h2o_ctm = 1;
            catch(copy_str(gas_optics->h2o_ctm_dir, h2o_ctm_dir, DIR_PATH_LEN));
        }
    }

    /*Prepare ozone continuum.*/
    gas_optics->use_o3_ctm = 0;
    if (o3_ctm_file != NULL)
    {
        if (strcmp(o3_ctm_file, "none") != 0)
        {
            gas_optics->use_o3_ctm = 1;
            catch(copy_str(gas_optics->o3_ctm_file, o3_ctm_file, DIR_PATH_LEN));
        }
    }

    /*Reserve memory.*/
    gmalloc(gas_optics->x, gas_optics->num_levels*NUM_MOLS, gas_optics->device);
    gmalloc(gas_optics->x_cfc, gas_optics->num_levels*NUM_CFCS, gas_optics->device);
    gmalloc(gas_optics->x_cia, gas_optics->num_levels*NUM_CIAS, gas_optics->device);
    gmalloc(gas_optics->n, gas_optics->num_layers, gas_optics->device);
    gmalloc(gas_optics->pavg, gas_optics->num_layers, gas_optics->device);
    gmalloc(gas_optics->tavg, gas_optics->num_layers, gas_optics->device);
    gmalloc(gas_optics->psavg, gas_optics->num_layers, gas_optics->device);
    gmalloc(gas_optics->ns, gas_optics->num_layers, gas_optics->device);
    gmalloc(gas_optics->linecenter, gas_optics->num_layers*MAX_NUM_LINES, gas_optics->device);
    gmalloc(gas_optics->snn, gas_optics->num_layers*MAX_NUM_LINES, gas_optics->device);
    gmalloc(gas_optics->gamma, gas_optics->num_layers*MAX_NUM_LINES, gas_optics->device);
    gmalloc(gas_optics->alpha, gas_optics->num_layers*MAX_NUM_LINES, gas_optics->device);
    if (gas_optics->device != HOST_ONLY)
    {
        gmalloc(gas_optics->p, gas_optics->num_levels, gas_optics->device);
        gmalloc(gas_optics->t, gas_optics->num_levels, gas_optics->device);
        gmalloc(gas_optics->tau, gas_optics->num_layers*gas_optics->grid.n,
                gas_optics->device);
    }

    /*Initialize TIPS.*/
    if (gas_optics->device != HOST_ONLY)
    {
        catch(inittips_d());
    }
    return GRTCODE_SUCCESS;
}


/*Free memory for the gas optics.*/
EXTERN int destroy_gas_optics(GasOptics_t * const gas_optics)
{
    not_null(gas_optics);
    int i;
    for (i=0; i<gas_optics->num_molecules; ++i)
    {
        catch(free_molecule(&(gas_optics->mols[i])));
    }
    for (i=0; i<gas_optics->num_cfcs; ++i)
    {
        catch(free_cfc_cross_sections(&(gas_optics->cfcs[i])));
    }
    for (i=0; i<gas_optics->num_cias; ++i)
    {
        catch(free_collision_induced_cross_sections(&(gas_optics->cia[i])));
    }
    catch(destroy_spectral_bins(&(gas_optics->bins)));
    gfree(gas_optics->x, gas_optics->device);
    gfree(gas_optics->x_cfc, gas_optics->device);
    gfree(gas_optics->x_cia, gas_optics->device);
    gfree(gas_optics->n, gas_optics->device);
    gfree(gas_optics->pavg, gas_optics->device);
    gfree(gas_optics->tavg, gas_optics->device);
    gfree(gas_optics->psavg, gas_optics->device);
    gfree(gas_optics->ns, gas_optics->device);
    gfree(gas_optics->linecenter, gas_optics->device);
    gfree(gas_optics->snn, gas_optics->device);
    gfree(gas_optics->gamma, gas_optics->device);
    gfree(gas_optics->alpha, gas_optics->device);
    if (gas_optics->device != HOST_ONLY)
    {
        gfree(gas_optics->p, gas_optics->device);
        gfree(gas_optics->t, gas_optics->device);
        gfree(gas_optics->tau, gas_optics->device);
    }
    int id;
    catch(molecule_hash(H2O, &id));
    if (gas_optics->use_h2o_ctm && is_active(gas_optics->molecule_bit_field, id))
    {
        catch(free_water_vapor_continuum_coefs(&(gas_optics->h2o_cc)));
    }
    catch(molecule_hash(O3, &id));
    if (gas_optics->use_o3_ctm && is_active(gas_optics->molecule_bit_field, id))
    {
        catch(free_ozone_continuum_coefs(&(gas_optics->o3_cc)));
    }
    return GRTCODE_SUCCESS;
}


/*Add a molecule.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int add_molecule(GasOptics_t * const gas_optics, int const molecule_id,
                        double const * const min_line_center,
                        double const * const max_line_center)
{
    not_null(gas_optics);
    int id;
    catch(molecule_hash(molecule_id, &id));
    if (is_active(gas_optics->molecule_bit_field, id))
    {
        char const *mesg = "molecule %d has already been added.";
        raise(GRTCODE_VALUE_ERR, mesg, molecule_id);
    }
    int index = gas_optics->num_molecules;
    (gas_optics->num_molecules)++;
    in_range(gas_optics->num_molecules, 1, NUM_MOLS);
    catch(activate(&(gas_optics->molecule_bit_field), id));
    double w0;
    if (min_line_center != NULL)
    {
        in_range(*min_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        w0 = *min_line_center;
    }
    else
    {
        w0 = gas_optics->grid.w0;
    }
    double wn;
    if (max_line_center != NULL)
    {
        in_range(*max_line_center, MIN_WAVENUMBER, MAX_WAVENUMBER);
        wn = *max_line_center;
    }
    else
    {
        wn = gas_optics->grid.wn;
    }
    min_check(wn, w0);
    catch(create_molecule(&(gas_optics->mols[index]), molecule_id, gas_optics->hitran_path,
                          w0, wn, gas_optics->num_layers, gas_optics->device));
    char const *mesg = "Using %s (%zu lines in range %e - %e [1/cm]).";
    log_mesg(mesg, gas_optics->mols[index].name, gas_optics->mols[index].line_params.num_lines,
             w0, wn);
    max_check(gas_optics->mols[index].line_params.num_lines, (uint64_t)(MAX_NUM_LINES-1));

    if (molecule_id == H2O && gas_optics->use_h2o_ctm)
    {
        /*Read in the water vapor continuum coefficients.*/
        mesg ="Using the %s continuum.";
        log_mesg(mesg, gas_optics->mols[index].name);
        catch(get_water_vapor_continuum_coefs(&(gas_optics->h2o_cc), gas_optics->h2o_ctm_dir,
                                              gas_optics->grid, gas_optics->device));
    }

    if (molecule_id == O3 && gas_optics->use_o3_ctm)
    {
        /*Read in the ozone continuum coefficients.*/
        mesg = "Using the %s continuum.";
        log_mesg(mesg, gas_optics->mols[index].name);
        catch(get_ozone_continuum_coefs(&(gas_optics->o3_cc), gas_optics->o3_ctm_file,
                                        gas_optics->grid, gas_optics->device));
    }
    return GRTCODE_SUCCESS;
}


/*Update a molecule's ppmv.*/
EXTERN int set_molecule_ppmv(GasOptics_t * const gas_optics, int const molecule_id,
                             fp_t const * const ppmv)
{
    not_null(gas_optics);
    not_null(ppmv);
    int index;
    catch(molecule_hash(molecule_id, &index));
    if (!is_active(gas_optics->molecule_bit_field, index))
    {
        char const *mesg = "molecule %d is not being used.";
        log_warn(mesg, molecule_id);
        return GRTCODE_SUCCESS;
    }
    fp_t a[gas_optics->num_levels];
    int i;
    for (i=0; i<gas_optics->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = index*gas_optics->num_levels;
    gmemcpy(&(gas_optics->x[offset]), a, gas_optics->num_levels, gas_optics->device,
            FROM_HOST);
    return GRTCODE_SUCCESS;
}


/*Add a CFC.*/
EXTERN int add_cfc(GasOptics_t * const gas_optics, int const cfc_id,
                   char const * const filepath)
{
    not_null(gas_optics);
    in_range(cfc_id, 0, NUM_CFCS);
    if (is_active(gas_optics->cfc_bit_field, cfc_id))
    {
        char const *mesg = "cfc %d has already been added.";
        raise(GRTCODE_VALUE_ERR, mesg, cfc_id);
    }
    int index = gas_optics->num_cfcs;
    (gas_optics->num_cfcs)++;
    in_range(gas_optics->num_cfcs, 1, NUM_CFCS);
    catch(activate(&(gas_optics->cfc_bit_field), cfc_id));

    /*Read in the CFC cross section values.*/
    catch(get_cfc_cross_sections(&(gas_optics->cfcs[index]), cfc_id, filepath,
                                 gas_optics->grid, gas_optics->device));
    char const *mesg = "Using CFC %s.";
    log_mesg(mesg, gas_optics->cfcs[index].name);
    return GRTCODE_SUCCESS;
}


/*Update a CFC's ppmv.*/
EXTERN int set_cfc_ppmv(GasOptics_t * const gas_optics, int const cfc_id,
                        fp_t const * const ppmv)
{
    not_null(gas_optics);
    not_null(ppmv);
    in_range(cfc_id, 0, NUM_CFCS);
    if (!is_active(gas_optics->cfc_bit_field, cfc_id))
    {
        char const *mesg = "CFC %d is not being used.";
        log_warn(mesg, cfc_id);
        return GRTCODE_SUCCESS;
    }
    fp_t a[gas_optics->num_levels];
    int i;
    for (i=0; i<gas_optics->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = cfc_id*gas_optics->num_levels;
    gmemcpy(&(gas_optics->x_cfc[offset]), a, gas_optics->num_levels, gas_optics->device,
            FROM_HOST);
    return GRTCODE_SUCCESS;
}


/*Activate collision-induced absorption between two species.*/
EXTERN int add_cia(GasOptics_t * const gas_optics, int const species1, int const species2,
                   char const * const filepath)
{
    not_null(gas_optics);
    int i;
    in_range(species1, 0, NUM_CIAS);
    in_range(species2, 0, NUM_CIAS);
    for (i=0; i<gas_optics->num_cias; ++i)
    {
        CollisionInducedAbsorption_t *m = &(gas_optics->cia[i]);
        if ((m->id[0] + m->id[1]) == (species1 + species2))
        {
            char const *mesg = "CIA with %s and %s is already active.";
            raise(GRTCODE_VALUE_ERR, mesg, m->name[0], m->name[1]);
        }
    }
    int index = gas_optics->num_cias;
    (gas_optics->num_cias)++;
    if (!is_active(gas_optics->cia_bit_field, species1))
    {
        catch(activate(&(gas_optics->cia_bit_field), species1));
    }
    if (!is_active(gas_optics->cia_bit_field, species2))
    {
        catch(activate(&(gas_optics->cia_bit_field), species2));
    }
    int id[2] = {species1, species2};
    catch(get_collision_induced_cross_sections(&(gas_optics->cia[index]), id,
                                               filepath, gas_optics->grid, gas_optics->device));
    char const *mesg = "Using collision-induced absorption between %s and %s.";
    log_info(mesg, gas_optics->cia[index].name[0], gas_optics->cia[index].name[1]);
    return GRTCODE_SUCCESS;
}


/*Update a CIA species' ppmv.*/
EXTERN int set_cia_ppmv(GasOptics_t * const gas_optics, int const cia_id,
                        fp_t const * const ppmv)
{
    not_null(gas_optics);
    not_null(ppmv);
    in_range(cia_id, 0, NUM_CIAS);
    if (!is_active(gas_optics->cia_bit_field, cia_id))
    {
        char const *mesg = "CIA %d is not being used.";
        log_warn(mesg, cia_id);
        return GRTCODE_SUCCESS;
    }
    fp_t a[MAX_NUM_LEVELS];
    int i;
    for (i=0; i<gas_optics->num_levels; ++i)
    {
        a[i] = ppmv[i]*1.e-6;
    }
    int offset = cia_id*gas_optics->num_levels;
    gmemcpy(&(gas_optics->x_cia[offset]), a, gas_optics->num_levels, gas_optics->device, FROM_HOST);
    return GRTCODE_SUCCESS;
}


/*Calcluate the total optical depth in each layer at each spectral grid point.*/
EXTERN int calculate_optical_depth(GasOptics_t * const gas_optics, fp_t * const pressure,
                                   fp_t * const temperature, Optics_t * const optics)
{
    not_null(gas_optics);
    not_null(pressure);
    not_null(temperature);
    not_null(optics);
    assert(gas_optics->device, optics->device);
    assert(gas_optics->num_layers, optics->num_layers);
    int same_grids;
    catch(compare_spectral_grids(&(gas_optics->grid), &(optics->grid), &same_grids));
    assert(same_grids, 1);
    fp_t const mbtoatm = 0.000986923f;
    fp_t p[MAX_NUM_LEVELS];
    int i;
    for (i=0; i<gas_optics->num_levels; ++i)
    {
        p[i] = pressure[i]*mbtoatm;
    }
    catch(launch(gas_optics, p, temperature, optics->tau));
    return GRTCODE_SUCCESS;
}


/*Get the number of molecules.*/
EXTERN int get_num_molecules(GasOptics_t const * const gas_optics, int * const n)
{
    not_null(gas_optics);
    not_null(n);
    *n = gas_optics->num_molecules;
    return GRTCODE_SUCCESS;
}
