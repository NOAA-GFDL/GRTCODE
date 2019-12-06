/** @file*/
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

#ifndef VERBOSITY_H_
#define VERBOSITY_H_


/** @brief Verbosity levels.*/
enum grtcode_verbosity
{
    GRTCODE_NONE, /**< Do not any messages.*/
    GRTCODE_ERROR, /**< Include error messages.*/
    GRTCODE_WARN, /**< Print warnings.*/
    GRTCODE_INFO /**< Print informational messages.*/
};


/** @brief Return a message for an input return code.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int grtcode_errstr(int const code, /**< Error code.*/
                          char * const buffer, /**< Buffer to hold error message.*/
                          int const buffer_size /**< Size of input buffer.*/
                         );


/** @brief Set the current verbosity level.*/
EXTERN void grtcode_set_verbosity(int const level /**< Verbosity level.*/
                                 );


/** @brief Retrieve the current verbosity level.
    @return The current verbosity level.*/
EXTERN int grtcode_verbosity();


#endif
