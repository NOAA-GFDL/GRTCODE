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

#ifndef CFCS_H_
#define CFCS_H_

#include <stdint.h>
#include "grtcode_utilities.h"


/** @brief Maximum size of CFC name.*/
#define CFC_NAME_LEN 16


/** @brief CFC identifiers.*/
typedef enum CfcId
{
    CFC11 = 0,
    CFC12,
    CFC113,
    CFC114,
    CFC115,
    HCFC22,
    HCFC141b,
    HCFC142b,
    HFC23,
    HFC125,
    HFC134a,
    HFC143a,
    HFC152a,
    HFC227ea,
    HFC245fa,
    CCl4,
    C2F6,
    CF4,
    CH2Cl2,
    NF3,
    SF6,
    NUM_CFCS
} CfcId_t;


/** @brief CFC cross-sections.*/
typedef struct CfcCrossSection
{
    fp_t *cross_section; /**< CFC cross-section [cm2].*/
    int id; /**< CFC id.*/
    char name[CFC_NAME_LEN]; /**< CFC name.*/
    uint64_t num_wpoints; /**< Spectral grid size.*/
    Device_t device; /**< Device id.*/
} CfcCrossSection_t;


#endif
