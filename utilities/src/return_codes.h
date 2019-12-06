/**  @file*/
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

#ifndef RETURN_CODES_H_
#define RETURN_CODES_H_


/** @brief Codes returned by GRTCODE functions.*/
enum grtcode_return_codes
{
    GRTCODE_SUCCESS = 0, /**< No error occurred.*/
    GRTCODE_INVALID_ERR, /**< Detected floating-point invalid.*/
    GRTCODE_DIVBYZERO_ERR, /**< Detected divided by zero.*/
    GRTCODE_OVERFLOW_ERR, /**< Detected overflow.*/
    GRTCODE_UNDERFLOW_ERR, /**< Detected underflow.*/
    GRTCODE_SENTINEL_ERR, /**< Entered unallowed code path.*/
    GRTCODE_NULL_ERR, /**< Attempted to dereference null pointer.*/
    GRTCODE_NON_NULL_ERR, /**< Attempted to malloc a non-null pointer.*/
    GRTCODE_RANGE_ERR, /**< Value is outside of allowed range.*/
    GRTCODE_VALUE_ERR, /**< Dectected illegal value.*/
    GRTCODE_COMPILER_ERR, /**< Requested feature not supported by current compiler.*/
    GRTCODE_IO_ERR, /**< Detected error while performing I/O.*/
    GRTCODE_GPU_ERR /**< Detected error on GPU device.*/
};


#endif
