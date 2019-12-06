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

#ifndef PARSE_CSV_H_
#define PARSE_CSV_H_


/** @brief Parse a csv file, assuming that the first line in the file contains
           headers for each of the columns.
    @return GRTCODE_SUCCESS or an error code.*/
int parse_csv(char const * const filepath, /**< csv file.*/
              int * const num_lines, /**< Number of lines in the file.*/
              int * const num_cols, /**< Number of columns in the file.*/
              int const ignore_headers, /**< Flag to ignore the first line in the file.*/
              char *** out /**< Array of values.*/
             );


#endif
