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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "cfcs.h"
#include "cfcs-internal.h"
#include "debug.h"
#include "grtcode_utilities.h"
#include "parse_csv.h"


/*Read in the cfc cross sections.*/
int get_cfc_cross_sections(CfcCrossSection_t *xsc, int const id, char const * const filepath,
                           SpectralGrid_t const grid, Device_t const device)
{
    not_null(xsc);
    not_null(filepath);

    /*Set CFC metadata.*/
    xsc->id = id;
    switch(id)
    {
        case CFC11:
            snprintf(xsc->name, CFC_NAME_LEN, "CFC-11");
            break;
        case CFC12:
            snprintf(xsc->name, CFC_NAME_LEN, "CFC-12");
            break;
        case CFC113:
            snprintf(xsc->name, CFC_NAME_LEN, "CFC-113");
            break;
        case CFC114:
            snprintf(xsc->name, CFC_NAME_LEN, "CFC-114");
            break;
        case CFC115:
            snprintf(xsc->name, CFC_NAME_LEN, "CFC-115");
            break;
        case HCFC22:
            snprintf(xsc->name, CFC_NAME_LEN, "HCFC-22");
            break;
        case HCFC141b:
            snprintf(xsc->name, CFC_NAME_LEN, "HCFC-141b");
            break;
        case HCFC142b:
            snprintf(xsc->name, CFC_NAME_LEN, "HCFC-142b");
            break;
        case HFC23:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-23");
            break;
        case HFC125:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-125");
            break;
        case HFC134a:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-134a");
            break;
        case HFC143a:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-143a");
            break;
        case HFC152a:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-152a");
            break;
        case HFC227ea:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-227ea");
            break;
        case HFC245fa:
            snprintf(xsc->name, CFC_NAME_LEN, "HFC-245fa");
            break;
        case CCl4:
            snprintf(xsc->name, CFC_NAME_LEN, "CCl4");
            break;
        case C2F6:
            snprintf(xsc->name, CFC_NAME_LEN, "C2F6");
            break;
        case CF4:
            snprintf(xsc->name, CFC_NAME_LEN, "CF4");
            break;
        case CH2Cl2:
            snprintf(xsc->name, CFC_NAME_LEN, "CH2Cl2");
            break;
        case NF3:
            snprintf(xsc->name, CFC_NAME_LEN, "NF3");
            break;
        case SF6:
            snprintf(xsc->name, CFC_NAME_LEN, "SF6");
            break;
        default:
            {char const *mesg = "unrecognized CFC id %d.";
            raise(GRTCODE_VALUE_ERR, mesg, id);}
    }

    /*Read in the data.*/
    int num_lines;
    int num_cols;
    char **buf;
    catch(parse_csv(filepath, &num_lines, &num_cols, 1, &buf));
    int const ncols_req = 2;
    if (num_cols != ncols_req)
    {
        char const *mesg = "The number of columns (%d) in file %s does not match"
                           " the expected number (%d).";
        raise(GRTCODE_VALUE_ERR, mesg, num_cols, filepath, ncols_req);
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
        xsc->cross_section = c;
    }
    else
    {
        gmalloc(xsc->cross_section, grid.n, device);
        gmemcpy(xsc->cross_section, c, grid.n, device, FROM_HOST);
        gfree(c, HOST_ONLY);
    }
    xsc->num_wpoints = grid.n;
    xsc->device = device;
    return GRTCODE_SUCCESS;
}


/*Free memory.*/
int free_cfc_cross_sections(CfcCrossSection_t *xsc)
{
    not_null(xsc);
    gfree(xsc->cross_section, xsc->device);
    return GRTCODE_SUCCESS;
}
