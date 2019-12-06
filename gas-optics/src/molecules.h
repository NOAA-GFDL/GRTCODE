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

#ifndef MOLECULES_H_
#define MOLECULES_H_

#include "grtcode_utilities.h"
#include "parse_HITRAN_file.h"


/** @brief Maximum size of molecule name.*/
#define MOL_NAME_LEN 8


/*Note: these numbers are assigned by the HITRAN database.*/
typedef enum HitranMoleculeId
{
    H2O = 1,
    CO2 = 2,
    O3 = 3,
    N2O = 4,
    CO = 5,
    CH4 = 6,
    O2 = 7,
    NO = 8,
    SO2 = 9,
    NO2 = 10,
    NH3 = 11,
    HNO3 = 12,
    OH = 13,
    HF = 14,
    HCl = 15,
    HBr = 16,
    HI = 17,
    ClO = 18,
    OCS = 19,
    H2CO = 20,
    HOCl = 21,
    N2 = 22,
    HCN = 23,
    CH3Cl = 24,
    H2O2 = 25,
    C2H2 = 26,
    C2H6 = 27,
    PH3 = 28,
    COF2 = 29,
    SF6_MOL = 30,
    H2S = 31,
    HCOOH = 32,
    HO2 = 33,
    O = 34,
    ClONO2 = 35,
    NOp = 36,
    HOBr = 37,
    C2H4 = 38,
    CH3OH = 39,
    CH3Br = 40,
    CH3CN = 41,
    CF4_MOL = 42,
    C4H2 = 43,
    HC3N = 44,
    H2 = 45,
    CS = 46,
    SO3 = 47,
    C2N2 = 48,
    COCl2 = 49,
    SO = 50,
    C3H4 = 51,
    CH3 = 52,
    CS2 = 53,
    NUM_MOLS = 53
} HitranMoleculeId_t;


/** @brief Molecule object.*/
typedef struct Molecule
{
    Device_t device; /**< Device id.*/
    int id; /**< HITRAN id.*/
    LineParams_t line_params; /**< Line parameters.*/
    fp_t mass; /**< Mass [g].*/
    char name[MOL_NAME_LEN]; /**< Name.*/
    int num_isotopologues; /**< Number of isotopologues.*/
    fp_t *q; /**< Total partition function (isotopologue, layer).*/
} Molecule_t;


#endif
