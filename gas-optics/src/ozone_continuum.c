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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "grtcode_utilities.h"
#include "ozone_continuum.h"
#include "ozone_continuum-internal.h"
#include "parse_csv.h"


/*Read in the ozone continuum coefficients.*/
int get_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc, char const * const o3_ctm_file,
                              SpectralGrid_t const grid, Device_t const device)
{
    not_null(cc);
    not_null(o3_ctm_file);

    /*Read in the data.*/
    char const *mesg = "Reading in ozone continuum coefficients from file %s.";
    log_info(mesg, o3_ctm_file);
    int num_lines;
    int num_cols;
    char **buf;
    catch(parse_csv(o3_ctm_file, &num_lines, &num_cols, 1, &buf));
    if (num_cols != 2)
    {
        mesg = "The number of columns (%d) in file %s does not match"
               " the expected number (%d).";
        raise(GRTCODE_VALUE_ERR, mesg, num_cols, o3_ctm_file, 2);
    }

    /*Convert the data from strings to floating point.*/
    fp_t *fbuf;
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

    /*Allocate space for the coefficient values at each wavenumber.*/
    fp_t *c;
    gmalloc(c, grid.n, HOST_ONLY);
    gmemset(c, 0, grid.n, HOST_ONLY);

    /*Interpolate to wavenumber grid.*/
    fp_t *x = &(fbuf[0]);
    fp_t *y = &(fbuf[num_lines]);
    catch(interpolate_to_grid(grid, x, y, (size_t)num_lines, c, linear_sample, NULL));
    gfree(fbuf, HOST_ONLY);
    if (device == HOST_ONLY)
    {
        cc->cross_section = c;
    }
    else
    {
        gmalloc(cc->cross_section, grid.n, device);
        gmemcpy(cc->cross_section, c, grid.n, device, FROM_HOST);
        gfree(c, HOST_ONLY);
    }
    cc->num_wpoints = grid.n;
    cc->device = device;
    return GRTCODE_SUCCESS;
}


/*Free memory used for the ozone continuum cross-sections.*/
int free_ozone_continuum_coefs(OzoneContinuumCoefs_t *cc)
{
    not_null(cc);
    gfree(cc->cross_section, cc->device);
    return GRTCODE_SUCCESS;
}
