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

#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "extern.h"
#include "gas_optics.h"
#include "optics.h"
#include "return_codes.h"
#include "solar_flux.h"
#include "spectral_grid.h"


enum StructTypes_t {
    GRID = 0,
    OPTICS,
    GAS_OPTICS,
    SOLAR_FLUX
};


/** @brief Malloc struct and associate input pointer.*/
/** @return GRTCODE_SUCCESS or an error code.*/
EXTERN int malloc_struct(void **p, int type)
{
    size_t s;
    switch (type)
    {
        case GRID:
            s = sizeof(SpectralGrid_t);
            break;
        case OPTICS:
            s = sizeof(Optics_t);
            break;
        case GAS_OPTICS:
            s = sizeof(GasOptics_t);
            break;
        case SOLAR_FLUX:
            s = sizeof(SolarFlux_t);
            break;
        default:
            {char *mesg = "unrecognized structure type %d.";
            raise(GRTCODE_VALUE_ERR, mesg, type);}
    }
    not_null(p);
    is_null(*p);
    *p = malloc(s);
    not_null(*p);
    return GRTCODE_SUCCESS;
}


/** @brief Free struct associated with input pointer.*/
/** @return GRTCODE_SUCCESS or an error code.*/
EXTERN int free_struct(void **p)
{
    gfree(*p, HOST_ONLY);
    *p = NULL;
    return GRTCODE_SUCCESS;
}


/** @brief Retrieve arrays from optics structure.*/
/** @return GRTCODE_SUCCESS or an error code.*/
EXTERN int optical_properties(Optics_t const * const optics, fp_t * const tau, fp_t * const omega,
                              fp_t * const g)
{
    not_null(optics);
    size_t n = (optics->num_layers)*(optics->grid.n);
    if (tau != NULL)
    {
        gmemcpy(tau, optics->tau, n, optics->device, FROM_DEVICE);
    }
    if (omega != NULL)
    {
        gmemcpy(omega, optics->omega, n, optics->device, FROM_DEVICE);
    }
    if (g != NULL)
    {
        gmemcpy(g, optics->g, n, optics->device, FROM_DEVICE);
    }
    return GRTCODE_SUCCESS;
}


/** @brief Get the spectral grid properties.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int spectral_grid_properties(SpectralGrid_t const * const grid, /**< Spectral grid.*/
                                    double * const w0, /**< Grid lower bound.*/
                                    uint64_t * const n, /**< Spectral grid size.*/
                                    double * const dw /**< Grid spacing.*/
                                   )
{
    not_null(grid);
    if (w0 != NULL)
    {
        *w0 = grid->w0;
    }
    if (n != NULL);
    {
        *n = grid->n;
    }
    if (dw != NULL)
    {
        *dw = grid->dw;
    }
    return GRTCODE_SUCCESS;
}


/** @brief Get the solar flux properties.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int solar_flux_properties(SolarFlux_t const * const solar, /**< Solar flux.*/
                                 fp_t * const flux /**< Flux.*/
                                )
{
    not_null(solar);
    not_null(flux);
    memcpy(flux, solar->incident_flux, sizeof(*(solar->incident_flux))*solar->n);
    return GRTCODE_SUCCESS;
}
