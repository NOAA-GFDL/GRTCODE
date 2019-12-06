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
#include <string.h>
#include "debug.h"
#include "extern.h"
#include "return_codes.h"
#include "verbosity.h"
#include "verbosity-internal.h"


static int verbosity = GRTCODE_NONE; /*Verbosity level.*/
static char error_buffer[4096]; /*Buffer containing current errors and backtrace.*/


/*Add an error to the error buffer.*/
EXTERN void append_to_error_buffer(char const * const mesg)
{
    char b[4096];
    snprintf(b, 4096, "%s", error_buffer);
    snprintf(error_buffer, 4096, "%s%s", b, mesg);
}


/*Copy the contents of the error buffer.*/
EXTERN void copy_error_buffer(char * const buffer, int const buffer_size)
{
    snprintf(buffer, buffer_size, "%s\n", error_buffer);
}


/*Return a message for an input return code.*/
EXTERN int grtcode_errstr(int const code, char * const buffer, int const buffer_size)
{
    not_null(buffer);
    min_check(buffer_size, 1);
    if (code == GRTCODE_SUCCESS)
    {
        snprintf(buffer, buffer_size, "No errors.");
    }
    else
    {
        copy_error_buffer(buffer, buffer_size);
    }
    return GRTCODE_SUCCESS;
}


/*Set the current verbosity level.*/
EXTERN void grtcode_set_verbosity(int const level)
{
    verbosity = level;
}


/*Retrieve the current verbosity level.*/
EXTERN int grtcode_verbosity()
{
    return verbosity;
}


/*Clear the error buffer.*/
EXTERN void reset_error_buffer()
{
    memset(error_buffer, '\0', 4096*sizeof(*error_buffer));
}
