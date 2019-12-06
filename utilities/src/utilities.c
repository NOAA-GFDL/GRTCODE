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

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "extern.h"
#include "floating_point_type.h"
#include "grtcode_config.h"
#include "return_codes.h"
#include "utilities.h"


/*Turn on bit in bit field.*/
int activate(uint64_t * const bit_field, int const index)
{
    not_null(bit_field);
    in_range(index, 0, 63);
    uint64_t const one = 1;
    *bit_field = (*bit_field) | (one << index);
    return GRTCODE_SUCCESS;
}


/*Calculate an angstrom exponent.*/
fp_t angstrom_exponent(fp_t tau1, fp_t tau2, fp_t lambda1, fp_t lambda2)
{
    fp_t const c = -1.;
    return c*LOG(tau1/tau2)/LOG(lambda1/lambda2);
}


/*Intepolate between two points using an "angstrom exponent".*/
EXTERN int angstrom_exponent_sample(fp_t const * const x, fp_t const * const y,
                                    fp_t const * const newx, fp_t * const newy, size_t n)
{
    size_t i;
    for (i=0; i<2; ++i)
    {
        if (y[i] <= 0.)
        {
            char const * mesg = "Cannot calculate the angstrom exponent because"
                                " y[%zu] <= 0 (%e)";
            raise(GRTCODE_VALUE_ERR, mesg, i, y[i]);
        }
    }
    fp_t const alpha = -1.*angstrom_exponent(y[1], y[0], x[0], x[1]);
    for (i=0; i<n; ++i)
    {
        newy[i] = y[0]*pow((x[0]/newx[i]), alpha);
    }
    return GRTCODE_SUCCESS;
}


/*Extrapolate using constant values.*/
int constant_extrapolation(fp_t const * const x, fp_t const * const y,
                           fp_t const * const newx, fp_t * const newy, size_t n)
{
    /*Note: the parameters x and new are unused, but are included so that this
      function matches the Sample1d_t typedef.  The following two lines exist
      only to quiet compiler warnings.*/
    (void)(x);
    (void)(newx);

    size_t i = 0;
    for (i=0; i<n; ++i)
    {
        newy[i] = y[0];
    }
    return GRTCODE_SUCCESS;
}


/*Copy a string into a buffer, checking its length.*/
int copy_str(char * const dest, char const * const src, size_t const len)
{
    if (strlen(src) > len)
    {
        char const *mesg = "input string (%s) is larger than the input buffer"
                           " and would be truncated.";
        raise(GRTCODE_VALUE_ERR, mesg, src);
    }
    snprintf(dest, len, "%s", src);
    return GRTCODE_SUCCESS;
}


/*Free malloced memory, with error checks.*/
EXTERN int free_ptr(void ** const p)
{
    not_null(p);
    not_null(*p);
    free(*p);
    *p = NULL;
    return GRTCODE_SUCCESS;
}


/*Integrate piece-wise between adjacent values in an array.*/
int integrate2(fp_t const * const x, fp_t const * const y, size_t n, fp_t * const s,
               Area1d_t area)
{
    not_null(x);
    not_null(y);
    not_null(s);
    not_null(area);
    if (n < 2)
    {
        char const * mesg = "at least two x points required (%zu given).";
        raise(GRTCODE_VALUE_ERR, mesg, n);
    }
    if (!monotonically_increasing(x, n))
    {
        char const * mesg = "x (%p) must be monotonically increasing.";
        raise(GRTCODE_VALUE_ERR, mesg, x);
    }
    *s = 0.;
    size_t i;
    for (i=0; i<(n-1); ++i)
    {
        *s += area(x+i, y+i);
    }
    return GRTCODE_SUCCESS;
}


/*Interpolate piece-wise between adjacent values in an array.*/
EXTERN int interpolate2(fp_t const * const x, fp_t const * const y, size_t n,
                        fp_t const * const newx, fp_t * const newy, size_t newn,
                        Sample1d_t interp, Sample1d_t extrap)
{
    not_null(x);
    not_null(y);
    not_null(newx);
    not_null(newy);
    if (n < 2)
    {
        char const * mesg = "at least two x points required (%zu given).";
        raise(GRTCODE_VALUE_ERR, mesg, n);
    }
    if (!monotonically_increasing(x, n))
    {
        char const * mesg = "x (%p) must be monotonically increasing.";
        raise(GRTCODE_VALUE_ERR, mesg, x);
    }
    if (!monotonically_increasing(newx, newn))
    {
        char const * mesg = "newx (%p) must be monotonically increasing.";
        raise(GRTCODE_VALUE_ERR, mesg, newx);
    }
    size_t i;
    for (i=0; i<newn; ++i)
    {
        if (newx[i] > x[0])
        {
            break;
        }
    }
    if (i > 0)
    {
        /*Handle all points where newx[:i-1] <= x[0].*/
        if (extrap != NULL)
        {
            catch(extrap(x, y, newx, newy, i));
        }
        if (i == newn)
        {
            /*All newx[:] <= x[0].*/
            return GRTCODE_SUCCESS;
        }
    }
    size_t j;
    for (j=0; j<n-1; ++j)
    {
        size_t k;
        for (k=i; k<newn; ++k)
        {
            if (newx[k] > x[j+1])
            {
                break;
            }
        }
        if (k > i)
        {
            /*Handle all points where x[j] < newx[i:k-1] <= x[j+1].*/
            catch(interp(&(x[j]), &(y[j]), &(newx[i]), &(newy[i]), k-i));
            i = k;
            if (i == newn)
            {
                /*All newx points have been handled.*/
                return GRTCODE_SUCCESS;
            }
        }
    }
    /*Handle all points where newx[i:] > x[n-1].*/
    if (extrap != NULL)
    {
        catch(extrap(&(x[n-2]), &(y[n-2]), &(newx[i]), &(newy[i]), newn-i));
    }
    return GRTCODE_SUCCESS;
}


/*Check if bit is turned on in bit field.*/
int is_active(uint64_t const bit_field, int const index)
{
    in_range(index, 0, (int)(CHAR_BIT*sizeof(bit_field)));
    uint64_t const one = 1;
    return bit_field & (one << index);
}


/*Linearly interpolate between two points.*/
EXTERN int linear_sample(fp_t const * const x, fp_t const * const y,
                         fp_t const * const newx, fp_t * const newy, size_t n)
{
    fp_t const m = (y[1] - y[0])/(x[1] - x[0]);
    fp_t const b = y[0] - m*x[0];
    size_t i;
    for (i=0; i<n; ++i)
    {
        newy[i] = m*newx[i] + b;
    }
    return GRTCODE_SUCCESS;
}


/*Malloc memory, with error checks.*/
EXTERN int malloc_ptr(void ** const p, size_t const num_bytes)
{
    not_null(p);
    *p = malloc(num_bytes);
    not_null(*p);
    return GRTCODE_SUCCESS;
}


/*Determine if the values in an array are monotonically increasing.*/
int monotonically_increasing(fp_t const * const x, size_t n)
{
  size_t i;
  for (i=0; i<(n-1); ++i)
  {
      if (x[i+1] <= x[i])
      {
          return 0;
      }
  }
  return 1;
}


/*Open a file, with error checks.*/
int open_file(FILE **file, char const * const name, char const * const mode)
{
    not_null(file);
    not_null(name);
    not_null(mode);
    *file = fopen(name, mode);
    if (*file == NULL)
    {
        char const *mesg = "failed to open file %s.";
        raise(GRTCODE_IO_ERR, mesg, name);
    }
    return GRTCODE_SUCCESS;
}


/*Convert a string to a double.*/
int to_double(char const * const s, double * const d)
{
    not_null(s);
    not_null(d);
    char *end;
    errno = 0;
    *d = strtod(s, &end);
    if ((errno == ERANGE && (*d == HUGE_VAL || *d == -1.*HUGE_VAL)) ||
        (*d == 0. && errno != 0))
    {
        char const *mesg = "the input string %s is out of range.";
        raise(GRTCODE_RANGE_ERR, mesg, s);
    }
    if (end == s)
    {
        char const *mesg = "invalid input string %s, expecting the string to"
                           " contain a floating point number.";
        raise(GRTCODE_VALUE_ERR, mesg, s);
    }
    return GRTCODE_SUCCESS;
}


/*Convert a double to fp_t.*/
int to_fp_t(double const d, fp_t * const f)
{
    not_null(f);
    if (sizeof(fp_t) == sizeof(double))
    {
        *f = d;
    }
    else if (sizeof(fp_t) == sizeof(float))
    {
        if (d >= -1.*FLT_MAX && d <= FLT_MAX)
        {
            *f = (fp_t)d;
        }
        else
        {
            char const *mesg = "input double value %le cannot be represented as a float.";
            raise(GRTCODE_RANGE_ERR, mesg, d);
        }
    }
    else
    {
        char const *mesg = "fp_t (size=%lu) must be represent either float or double.";
        raise(GRTCODE_VALUE_ERR, mesg, sizeof(fp_t));
    }
    return GRTCODE_SUCCESS;
}


/*Convert a string to an integer.*/
int to_int(char const * const s, int * const i)
{
    not_null(s);
    not_null(i);
    char *end;
    errno = 0;
    long n = strtol(s, &end, 10);
    if ((errno == ERANGE && (n == LONG_MAX || n == LONG_MIN)) ||
        (n == 0 && errno != 0))
    {
        char const *mesg = "the input string %s is out of range.";
        raise(GRTCODE_RANGE_ERR, mesg, s);
    }
    if (end == s || errno == EINVAL)
    {
        char const *mesg = "invalid input string %s, expecting the string to"
                           " contain an integer.";
        raise(GRTCODE_VALUE_ERR, mesg, s);
    }
    if (n >= INT_MIN && n <= INT_MAX)
    {
        *i = n;
    }
    else
    {
        char const *mesg = "input string %s cannot be represented as an int.";
        raise(GRTCODE_VALUE_ERR, mesg, s);
    }
    return GRTCODE_SUCCESS;
}


/*Calculate the integral between two points using the trapezoid rule.*/
fp_t trapezoid(fp_t const * const x, fp_t const * const y)
{
    fp_t const half = 0.5;
    return half*(y[0] + y[1])*(x[1] - x[0]);
}
