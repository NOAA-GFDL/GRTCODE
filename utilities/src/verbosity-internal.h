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

#ifndef VERBOSITY_INTERNAL_H_
#define VERBOSITY_INTERNAL_H_

#include "extern.h"


/** @brief Add an error to the error buffer.*/
EXTERN void append_to_error_buffer(char const * const mesg /**< Error message.*/
                                  );


/** @brief Copy the contents of the error buffer.*/
EXTERN void copy_error_buffer(char * const buffer, /**< Input buffer.*/
                              int const buffer_size /**< Size of the input buffer.*/
                             );


/** @brief Clear the error buffer.*/
EXTERN void reset_error_buffer();


#endif
