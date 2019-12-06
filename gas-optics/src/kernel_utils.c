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
#include "kernel_utils.h"


/*Get the indices of an array that bracket a value.  The input array must be sorted.*/
HOST DEVICE int bracket(uint64_t const array_size, fp_t const * const array, fp_t const val,
                        uint64_t * const left, uint64_t * const right)
{
    not_null(array);
    not_null(left);
    not_null(right);
    min_check(array_size, 1);
    uint64_t l = 0;
    uint64_t r = array_size - 1;
    if (val < array[l] || val > array[r])
    {
        *left = l;
        *right = r;
        return GRTCODE_RANGE_ERR;
    }
    if (array[l] == val)
    {
        r = l;
    }
    else if (array[r] == val)
    {
        l = r;
    }
    else
    {
        while (r - l > 1)
        {
            uint64_t mid = l + (r - l)/2;
            if (array[mid] == val)
            {
                l = mid;
                r = mid;
                break;
            }
            else if (val > array[mid])
            {
                l = mid;
            }
            else
            {
                r = mid;
            }
            if (l > r || l == r)
            {
                sentinel();
            }
        }
    }
    *left = l;
    *right = r;
    return GRTCODE_SUCCESS;
}


/*Perform a quadratic interpolation, but set all values less than zero equal to zero.*/
HOST DEVICE int bin_quad_interp(fp_t const * const x, fp_t const * const y,
                                uint64_t const left, uint64_t const right, fp_t const w0,
                                fp_t const dw, fp_t * const tau)
{
    not_null(x);
    not_null(y);
    not_null(tau);
    uint64_t j;
    for (j=left; j<=right; ++j)
    {
        fp_t w = w0 + j*dw;
        fp_t t = (w - x[1])*(w - x[2])*y[0]/((x[0] - x[1])*(x[0] - x[2])) +
                 (w - x[0])*(w - x[2])*y[1]/((x[1] - x[0])*(x[1] - x[2])) +
                 (w - x[0])*(w - x[1])*y[2]/((x[2] - x[0])*(x[2] - x[1]));
        if (t < 0.f)
        {
            t = 0.f;
        }
        tau[j] += t;
    }
    return GRTCODE_SUCCESS;
}


/*Copy optical depths from the coarse mesh to the fine mesh in a bin.*/
HOST DEVICE int bin_no_interp(uint64_t const left, uint64_t const right,
                              fp_t const * const taub, fp_t * const tau)
{
    not_null(taub);
    not_null(tau);
    uint64_t j;
    for (j=left; j<=right; ++j)
    {
        tau[j] += taub[j - left];
    }
    return GRTCODE_SUCCESS;
}
