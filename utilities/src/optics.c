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

#include <stdint.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "grtcode_config.h"
#include "optics.h"
#include "return_codes.h"
#include "spectral_grid.h"


/*Add optics objects together.*/
static int add_optics_objects(uint64_t const n, /*Size of output g, omega, and tau arrays.*/
                              int const num_optics, /*Number of mechanisms.*/
                              fp_t const * const g_in, /*Input asymmetric factor (n, mechanism).*/
                              fp_t const * const omega_in, /*Input single-scatter albedo (n, mechanism).*/
                              fp_t const * const tau_in, /*Input optical depth (n, mechanism).*/
                              fp_t * const g_out, /*Output asymmetric factor (n).*/
                              fp_t * const omega_out, /*Output single-scatter albedo (n).*/
                              fp_t * const tau_out /*Output optical depth (n).*/
                             );


#ifdef __NVCC__
/*Add optics objects together.*/
__global__ static void add_optics_objects_d(uint64_t const n, /*Size of output g, omega, and tau arrays.*/
                                            int const num_optics, /*Number of mechanisms.*/
                                            fp_t const * const g_in, /*Input asymmetric factor (n, mechanism).*/
                                            fp_t const * const omega_in, /*Input single-scatter albedo (n, mechanism).*/
                                            fp_t const * const tau_in, /*Input optical depth (n, mechanism).*/
                                            fp_t * const g_out, /*Output asymmetric factor (n).*/
                                            fp_t * const omega_out, /*Output single-scatter albedo (n).*/
                                            fp_t * const tau_out /*Output optical depth (n).*/
                                           );
#endif


/*Sample optical properties from another optics object.*/
static int sample_optics_object(uint64_t const n, /*Output spectral grid size.*/
                                uint64_t const factor, /*Input spectral grid size divided
                                                         by the output spectral grid size.*/
                                fp_t * const tau, /*Optical depth (layer, n).*/
                                fp_t * const omega, /*Single-scatter albedo (layer, n).*/
                                fp_t * const g, /*Asymmetric factor (layer, n).*/
                                fp_t const * const tau_in, /*Optical depth (layer, n*factor).*/
                                fp_t * const omega_in, /*Single-scatter albedo (layer, n*factor).*/
                                fp_t const * const g_in /*Asymmetric factor (layer, n*factor).*/
                               );


#ifdef __NVCC__
/*Sample optical properties from another optics object.*/
__global__ static void sample_optics_object_d(uint64_t const n, /*Output spectral grid size.*/
                                              uint64_t const factor, /*Input spectral grid size divided
                                                                       by the output spectral grid size.*/
                                              fp_t * const tau, /*Optical depth (layer, n).*/
                                              fp_t * const omega, /*Single-scatter albedo (layer, n).*/
                                              fp_t * const g, /*Asymmetric factor (layer, n).*/
                                              fp_t const * const tau_in, /*Optical depth (layer, n*factor).*/
                                              fp_t * const omega_in, /*Single-scatter albedo (layer, n*factor).*/
                                              fp_t const * const g_in /*Asymmetric factor (layer, n*factor).*/
                                             );
#endif


/*Add optical properties together.*/
EXTERN int add_optics(Optics_t const * const * const optics, int const num_optics,
                      Optics_t * const result)
{
    not_null(optics);
    not_null(result);
    min_check(num_optics, 1);
    Optics_t const *o = optics[0];
    not_null(o);
    int j;
    for (j=1; j<num_optics; ++j)
    {
        not_null(optics[j]);
        int ok;
        catch(optics_compatible(optics[j], o, &ok));
        if (ok == 0)
        {
            char const *mesg = "input optics objects (%p, %p) are incompatible.";
            raise(GRTCODE_VALUE_ERR, mesg, o, optics[j]);
        }
    }
    catch(create_optics(result, o->num_layers, &(o->grid), &(o->device)));
    fp_t *g = NULL;
    fp_t *omega = NULL;
    fp_t *tau = NULL;
    uint64_t n = o->num_layers*o->grid.n;
    gmalloc(g, n*num_optics, o->device);
    gmalloc(omega, n*num_optics, o->device);
    gmalloc(tau, n*num_optics, o->device);
    for (j=0; j<num_optics; ++j)
    {
        gmemcpy(&(g[j*n]), optics[j]->g, n, o->device, FROM_HOST);
        gmemcpy(&(omega[j*n]), optics[j]->omega, n, o->device, FROM_HOST);
        gmemcpy(&(tau[j*n]), optics[j]->tau, n, o->device, FROM_HOST);
    }
    glaunch(add_optics_objects, n, o->device, n, num_optics, g, omega, tau, result->g,
            result->omega, result->tau);
    gfree(g, o->device);
    gfree(omega, o->device);
    gfree(tau, o->device);
    return GRTCODE_SUCCESS;
}


/*Add optics objects together.*/
static int add_optics_objects(uint64_t const n, int const num_optics, fp_t const * const g_in,
                              fp_t const * const omega_in, fp_t const * const tau_in,
                              fp_t * const g_out, fp_t * const omega_out, fp_t * const tau_out)
{
    uint64_t i;
#pragma omp parallel for default(shared) private(i)
    for (i=0; i<n; ++i)
    {
        int j;
        for (j=0; j<num_optics; ++j)
        {
            uint64_t offset = j*n + i;
            g_out[i] += g_in[offset]*omega_in[offset]*tau_in[offset];
            omega_out[i] += omega_in[offset]*tau_in[offset];
            tau_out[i] += tau_in[offset];
        }
        g_out[i] /= omega_out[i];
        omega_out[i] /= tau_out[i];
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/*Add optics objects together.*/
__global__ static void add_optics_objects_d(uint64_t const n, int const num_optics,
                                            fp_t const * const g_in, fp_t const * const omega_in,
                                            fp_t const * const tau_in, fp_t * const g_out,
                                            fp_t * const omega_out, fp_t * const tau_out)
{
    uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n)
    {
        int j;
        for (j=0; j<num_optics; ++j)
        {
            uint64_t offset = j*n + i;
            g_out[i] += g_in[offset]*omega_in[offset]*tau_in[offset];
            omega_out[i] += omega_in[offset]*tau_in[offset];
            tau_out[i] += tau_in[offset];
        }
        g_out[i] /= omega_out[i];
        omega_out[i] /= tau_out[i];
    }
    return;
}
#endif


/*Reserve memory for the optics.*/
EXTERN int create_optics(Optics_t * const optics, int const num_layers, 
                         SpectralGrid_t const * const grid, Device_t const * const device)
{
    not_null(optics);
    not_null(grid);
    not_null(device);
    in_range(num_layers, MIN_NUM_LAYERS, MAX_NUM_LAYERS);
    optics->num_layers = num_layers;
    char const *mesg = "Initializing optics object:\nDevice: %d"
                       "\nAtmospheric column properties:\n\t"
                       "number of levels: %d\n\tnumber of layers: %d";
    log_info(mesg, *device, num_layers+1, num_layers);
    optics->grid = *grid;
    optics->device = *device;
    uint64_t n = num_layers*(optics->grid.n);
    gmalloc(optics->g, n, optics->device);
    gmemset(optics->g, 0, n, optics->device);
    gmalloc(optics->omega, n, optics->device);
    gmemset(optics->omega, 0, n, optics->device);
    gmalloc(optics->tau, n, optics->device);
    gmemset(optics->tau, 0, n, optics->device);
    return GRTCODE_SUCCESS;
}


/*Free memory for the optics.*/
EXTERN int destroy_optics(Optics_t * const optics)
{
    not_null(optics);
    gfree(optics->g, optics->device);
    gfree(optics->omega, optics->device);
    gfree(optics->tau, optics->device);
    return GRTCODE_SUCCESS;
}


/*Determine if two optics objects are compatible.*/
EXTERN int optics_compatible(Optics_t const * const one, Optics_t const * const two,
                             int * const result)
{
    not_null(one);
    not_null(two);
    not_null(result);
    int same_grids;
    catch(compare_spectral_grids(&(one->grid), &(two->grid), &same_grids));
    if ((one->num_layers == two->num_layers) && (same_grids == 1) &&
        (one->device == two->device))
    {
        *result = 1;
    }
    else
    {
        *result = 0;
    }
    return GRTCODE_SUCCESS;
}


/*Sample optical properties from another optics object.*/
EXTERN int sample_optics(Optics_t * const dest, Optics_t const * const source,
                         double const * const w0, double const * const wn)
{
    not_null(dest);
    not_null(source);
    if (dest->device != source->device)
    {
        char const *mesg = "Device ids (%d, %d) must match.";
        raise(GRTCODE_VALUE_ERR, mesg, dest->device, source->device);
    }
    if (dest->num_layers != source->num_layers)
    {
        char const *mesg = "Number of layers (%d, %d) must match.";
        raise(GRTCODE_VALUE_ERR, mesg, dest->num_layers, source->num_layers);
    }
    fp_t lower_bound;
    uint64_t lower_bound_index_d;
    uint64_t lower_bound_index_s;
    if (w0 == NULL)
    {
        lower_bound = dest->grid.w0;
        lower_bound_index_d = 0;
    }
    else
    {
        catch(grid_point_index(dest->grid, *w0, &lower_bound_index_d));
        lower_bound = *w0;
    }
    catch(grid_point_index(source->grid, lower_bound, &lower_bound_index_s));
    fp_t upper_bound;
    uint64_t upper_bound_index_d;
    uint64_t upper_bound_index_s;
    if (wn == NULL)
    {
        upper_bound = dest->grid.wn;
        upper_bound_index_d = dest->grid.n;
    }
    else
    {
        catch(grid_point_index(dest->grid, *wn, &upper_bound_index_d));
        upper_bound = *wn;
    }
    catch(grid_point_index(source->grid, upper_bound, &upper_bound_index_s));
    min_check(upper_bound, lower_bound);
    uint64_t const n_d = upper_bound_index_d - lower_bound_index_d + 1;
    uint64_t const n_s = upper_bound_index_s - lower_bound_index_s + 1;
    if ((n_d > n_s) || ((n_s - 1) % (n_d - 1) != 0))
    {
        char const *mesg = "New grid must be a subdomain of the original (%p).";
        raise(GRTCODE_VALUE_ERR, mesg, &(source->grid));
    }
    uint64_t const factor = (n_s - 1)/(n_d - 1);
    int i;
    for (i=0; i<dest->num_layers; ++i)
    {
        fp_t *tau_d = dest->tau + i*dest->grid.n + lower_bound_index_d;
        fp_t *omega_d = dest->omega + i*dest->grid.n + lower_bound_index_d;
        fp_t *g_d = dest->g + i*dest->grid.n + lower_bound_index_d;
        fp_t *tau_s = source->tau + i*source->grid.n + lower_bound_index_s;
        fp_t *omega_s = source->omega + i*source->grid.n + lower_bound_index_s;
        fp_t *g_s = source->g + i*source->grid.n + lower_bound_index_s;
        glaunch(sample_optics_object, n_d, dest->device, n_d, factor, tau_d, omega_d,
                g_d, tau_s, omega_s, g_s);
    }
    return GRTCODE_SUCCESS;
}


/*Sample optical properties from another optics object.*/
static int sample_optics_object(uint64_t const n, uint64_t const factor,
                                fp_t * const tau, fp_t * const omega, fp_t * const g,
                                fp_t const * const tau_in, fp_t * const omega_in,
                                fp_t const * const g_in)
{
    uint64_t i;
#pragma omp parallel for default(shared) private(i)
    for (i=0; i<n; ++i)
    {
        uint64_t offset = i*factor;
        tau[i] = tau_in[offset];
        omega[i] = omega_in[offset];
        g[i] = g_in[offset];
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/*Sample optical properties from another optics object.*/
__global__ static void sample_optics_object_d(uint64_t const n, uint64_t const factor,
                                              fp_t * const tau, fp_t * const omega,
                                              fp_t * const g, fp_t const * const tau_in,
                                              fp_t * const omega_in, fp_t const * const g_in)
{
    uint64_t i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n)
    {
        uint64_t offset = i*factor;
        tau[i] = tau_in[offset];
        omega[i] = omega_in[offset];
        g[i] = g_in[offset];
    }
    return;
}
#endif


/*Update optical properties.*/
EXTERN int update_optics(Optics_t * const optics, fp_t const * const tau,
                         fp_t const * const omega, fp_t const * const g)
{
    not_null(optics);
    not_null(tau);
    not_null(omega);
    not_null(g);
    uint64_t n = optics->num_layers*optics->grid.n;
    gmemcpy(optics->tau, tau, n, optics->device, FROM_HOST);
    gmemcpy(optics->omega, omega, n, optics->device, FROM_HOST);
    gmemcpy(optics->g, g, n, optics->device, FROM_HOST);
    return GRTCODE_SUCCESS;
}
