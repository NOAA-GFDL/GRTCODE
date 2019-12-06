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
#include "collision_induced_absorption.h"
#include "collision_induced_absorption-internal.h"
#include "debug.h"
#include "grtcode_utilities.h"
#include "parse_csv.h"


/*Read in the collision-induced absorption cross sections.*/
int get_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia,
                                         int const id[2], char const * const filepath,
                                         SpectralGrid_t const grid, Device_t const device)
{
    not_null(cia);
    not_null(filepath);

    /*Set CIA metadata.*/
    uint64_t j;
    for (j=0; j<2; ++j)
    {
        uint64_t offset = j*CIA_NAME_LEN;
        switch(id[j])
        {
            case CIA_N2:
                snprintf(&(cia->name_buf[offset]), CIA_NAME_LEN, "N2");
                break;
            case CIA_O2:
                snprintf(&(cia->name_buf[offset]), CIA_NAME_LEN, "O2");
                break;
            default:
                {char const *mesg = "unrecognized CIA id %d.";
                raise(GRTCODE_VALUE_ERR, mesg, id[j]);}
        }
        cia->id[j] = id[j];
        cia->name[j] = &(cia->name_buf[offset]);
    }

    /*Read in the data.*/
    int num_vals = 1;
    char const *mesg = "Reading in collision-induced absorption cross sections from file %s.";
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
    for (j=0; j<data_size; ++j)
    {
        double d;
        catch(to_double(buf[j], &d));
        catch(to_fp_t(d, &(fbuf[j])));
        gfree(buf[j], HOST_ONLY);
    }
    gfree(buf, HOST_ONLY);

    /*Allocate space for the cross section values at each wavenumber.*/
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
        cia->cross_section = c;
    }
    else
    {
        gmalloc(cia->cross_section, grid.n, device);
        gmemcpy(cia->cross_section, c, grid.n, device, FROM_HOST);
        gfree(c, HOST_ONLY);
    }
    cia->num_wpoints = grid.n;
    cia->device = device;
    return GRTCODE_SUCCESS;
}


/*Free memory reserved for collision-induced absorption cross sections.*/
int free_collision_induced_cross_sections(CollisionInducedAbsorption_t * const cia)
{
    not_null(cia);
    gfree(cia->cross_section, cia->device);
    return GRTCODE_SUCCESS;
}
