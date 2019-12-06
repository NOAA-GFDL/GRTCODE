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

#ifndef PARSE_HITRAN_FILE_H_
#define PARSE_HITRAN_FILE_H_

#include <stdint.h>
#include "grtcode_utilities.h"


/** @brief Molecular line parameters.*/
typedef struct LineParams
{
    fp_t *d; /**< Foreign-broadened pressure shift [cm-1 atm-1] of the transition
                  frequency at 296 K and 1 atm (line).*/
    Device_t device; /**< Device id.*/
    fp_t *en; /**< Lower state energy of the transition [cm-1] (line).*/
    int *iso; /**< Isotopologue id (line).*/
    fp_t *n; /**< Coefficient of temperature dependence of the foreign-broadened halfwidth
                  at half maximum (line).*/
    uint64_t num_lines; /**< Number of lines.*/
    fp_t *snn; /**< Adjusted line strength [cm-1] at 296K (line).*/
    fp_t *vnn; /**< Line transition frequency [cm-1] (line).*/
    fp_t *yair; /**< Foreign-broadened halfwidth at half maximum [cm-1 atm-1] at 296 K and
                     1 atm (line).*/
    fp_t *yself; /**< Self-broadened halfwidth at half maximum [cm-1 amt-1] at 296 K and
                      1 atm (line).*/
} LineParams_t;


#endif
