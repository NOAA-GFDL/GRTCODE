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
#include "grtcode_utilities.h"
#include "parse_csv.h"
#include "solar_flux.h"


/*Read in data for the solar flux.*/
EXTERN int create_solar_flux(SolarFlux_t * const solar_flux, SpectralGrid_t const * const grid,
                             char const * const filepath)
{
    not_null(solar_flux);
    not_null(grid);
    not_null(filepath);

    /*Read in the data.*/
    char *mesg = "Reading in solar flux values from file %s.";
    log_info(mesg, filepath);
    int num_vals = 1;
    int num_lines;
    int num_cols;
    char **buf;
    catch(parse_csv(filepath, &num_lines, &num_cols, 1, &buf));
    if ((num_vals + 1) != num_cols)
    {
        mesg = "The number of columns (%d) in file %s does not match"
               " the expected number (%d).";
        raise(GRTCODE_VALUE_ERR, mesg, num_cols, filepath, num_vals+1);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf = NULL;
    uint64_t data_size = num_lines*num_cols;
    gmalloc(fbuf, data_size, HOST_ONLY);
    uint64_t j;
    for (j=0; j<data_size; ++j)
    {
        double d;
        catch(to_double(buf[j], &d));
        catch(to_fp_t(d, &(fbuf[j])));
        gfree(buf[j], HOST_ONLY);
    }
    gfree(buf, HOST_ONLY);

    /*Allocate space for the solar fluxes.*/
    solar_flux->grid = *grid;
    fp_t *c = NULL;
    gmalloc(c, grid->n, HOST_ONLY);
    gmemset(c, 0, grid->n, HOST_ONLY);
    fp_t *w;
    catch(grid_points(*grid, &w, HOST_ONLY));

    /*Interpolate to spectral grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    catch(interpolate2(x, y, (size_t)num_lines, w, c, (size_t)grid->n, linear_sample, NULL));
    gfree(fbuf, HOST_ONLY);

    /*Integrate fluxes over spectral grid.*/
    fp_t total_flux;
    catch(integrate2(w, c, grid->n, &total_flux, trapezoid));
    gfree(w, HOST_ONLY);

    /*Adjust the spectral.*/
    for (j=0; j<grid->n; ++j)
    {
        c[j] /= total_flux;
    }
    solar_flux->incident_flux = c;
    solar_flux->n = grid->n;
    return GRTCODE_SUCCESS;
}


/*Free memory for the solar flux.*/
int destroy_solar_flux(SolarFlux_t * const solar_flux)
{
    not_null(solar_flux);
    gfree(solar_flux->incident_flux, HOST_ONLY);
    return GRTCODE_SUCCESS;
}
