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
#include "debug.h"
#include "device.h"
#include "extern.h"
#include "floating_point_type.h"
#include "grtcode_config.h"
#include "return_codes.h"
#include "spectral_grid.h"
#include "utilities.h"


/*Determine if two spectral grids are the same.*/
EXTERN int compare_spectral_grids(SpectralGrid_t const * const one,
                                  SpectralGrid_t const * const two, int * const result)
{
    not_null(one);
    not_null(two);
    not_null(result);
    if ((one->w0 == two->w0) && (one->wn == two->wn) && (one->dw == two->dw))
    {
        *result = 1;
    }
    else
    {
        *result = 0;
    }
    return GRTCODE_SUCCESS;
}


/*Initialize a spectral grid.*/
EXTERN int create_spectral_grid(SpectralGrid_t * const grid, double const w0, double const wn,
                                double const dw)
{
    not_null(grid);
    in_range(w0, MIN_WAVENUMBER, MAX_WAVENUMBER);
    grid->w0 = w0;
    in_range(wn, w0 + epsilon_, MAX_WAVENUMBER);
    grid->wn = wn;
    in_range(dw, MIN_RESOLUTION, MAX_RESOLUTION);
    grid->dw = dw;
    grid->n = ceil((wn-w0)/dw) + 1.;
    char const *mesg = "Spectral grid properties:\n\tlower bound: %e [1/cm]\n\t"
                       "upper bound: %e [1/cm]\n\tresolution: %e [1/cm]\n\t"
                       "total size: %zu grid points";
    log_info(mesg, w0, wn, dw, grid->n);
    return GRTCODE_SUCCESS;
}


/*Get the index of a grid point.*/
EXTERN int grid_point_index(SpectralGrid_t const grid, double const w, uint64_t * const index)
{
    not_null(index);
    in_range(w, grid.w0, grid.wn);
    *index = (uint64_t)(round((w - grid.w0)/grid.dw));
    double tolerance = grid.dw*1.e-5;
    if (fabs(grid.w0 + (*index)*grid.dw - w) > tolerance)
    {
        char const *mesg = "value %e not located on grid (%p).";
        raise(GRTCODE_VALUE_ERR, mesg, w, &grid);
    }
    return GRTCODE_SUCCESS;
}


/*Create an array of all grid points.*/
int grid_points(SpectralGrid_t const grid, fp_t **buffer, Device_t const device)
{
    fp_t *w;
    gmalloc(w, grid.n, device);
    uint64_t i;
    for (i=0; i<grid.n; ++i)
    {
        w[i] = grid.w0 + i*grid.dw;
    }
    *buffer = w;
    return GRTCODE_SUCCESS;
}


/*Interpolate an array onto the input spectral grid.*/
EXTERN int interpolate_to_grid(SpectralGrid_t const grid, fp_t const * const x,
                               fp_t const * const y, size_t const n,
                               fp_t * const newy, Sample1d_t interp,
                               Sample1d_t extrap)
{
    fp_t *w;
    catch(grid_points(grid, &w, HOST_ONLY));
    catch(interpolate2(x, y, n, w, newy, (size_t)grid.n, interp, extrap));
    gfree(w, HOST_ONLY);
    return GRTCODE_SUCCESS;
}
