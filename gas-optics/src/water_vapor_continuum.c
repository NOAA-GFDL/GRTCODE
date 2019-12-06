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
#include "debug.h"
#include "grtcode_utilities.h"
#include "parse_csv.h"
#include "water_vapor_continuum.h"
#include "water_vapor_continuum-internal.h"


/*Read in the water vapor continuum coefficients.*/
int get_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc,
                                    char const * const h2o_ctm_dir,
                                    SpectralGrid_t const grid,
                                    Device_t const device)
{
    not_null(cc);
    not_null(h2o_ctm_dir);
    char *filepath = NULL;
    size_t s = strlen(h2o_ctm_dir) + 64;
    gmalloc(filepath, s, HOST_ONLY);
    gmalloc(cc->coefs, NUM_COEFS, HOST_ONLY);
    int i;
    for (i=0; i<NUM_COEFS; ++i)
    {
        int num_vals;
        switch (i)
        {
            case MTCKD25_F296:
                snprintf(filepath, s, "%s/296MTCKD25_F.csv", h2o_ctm_dir);
                num_vals = 1;
                break;
            case MTCKD25_S296:
                snprintf(filepath, s, "%s/296MTCKD25_S.csv", h2o_ctm_dir);
                num_vals = 1;
                break;
            case CKDF:
                snprintf(filepath, s, "%s/CKDF.csv", h2o_ctm_dir);
                num_vals = 3;
                break;
            case CKDS:
                snprintf(filepath, s, "%s/CKDS.csv", h2o_ctm_dir);
                num_vals = 3;
                break;
            default:
                sentinel();
        }

        /*Read in the data.*/
        char const *mesg = "Reading in water vapor continuum coefficients from file %s.";
        log_info(mesg, filepath);
        int num_lines;
        int num_cols;
        char **buf;
        catch(parse_csv(filepath, &num_lines, &num_cols, 1, &buf));
        if ((num_vals+1) != num_cols)
        {
            mesg = "The number of columns (%d) in file %s does not match"
                   " the expected number (%d).";
            raise(GRTCODE_VALUE_ERR, mesg, num_cols, filepath, num_vals+1);
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
            cc->coefs[i] = c;
        }
        else
        {
            gmalloc(cc->coefs[i], grid.n, device);
            gmemcpy(cc->coefs[i], c, grid.n, device, FROM_HOST);
            gfree(c, HOST_ONLY);
        }
    }
    cc->num_wpoints = grid.n;
    cc->device = device;
    gfree(filepath, HOST_ONLY);
    return GRTCODE_SUCCESS;
}


/*Free memory used for the continuum coefficients.*/
int free_water_vapor_continuum_coefs(WaterVaporContinuumCoefs_t *cc)
{
    not_null(cc);
    int i;
    for (i=0; i<NUM_COEFS; ++i)
    {
        gfree(cc->coefs[i], cc->device);
    }
    gfree(cc->coefs, HOST_ONLY);
    return GRTCODE_SUCCESS;
}
