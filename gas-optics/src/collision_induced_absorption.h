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

#ifndef COLLISION_INDUCED_ABSORPTION_H_
#define COLLISION_INDUCED_ABSORPTION_H_

#include <stdint.h>
#include "grtcode_utilities.h"


/** @brief Maximum size of collision-induced absorption species name.*/
#define CIA_NAME_LEN 8


/*Number of possible CIA combinations (i.e. N2-N2, N2-O2, O2-O2).*/
#define MAX_NUM_CIAS 3


/*Collision-induced absorption identifiers.*/
typedef enum CiaId
{
    CIA_N2 = 0,
    CIA_O2,
    NUM_CIAS,
} CiaId_t;


/** @brief Collision-induced absorption object.*/
typedef struct CollisionInducedAbsorption
{
    int id[2]; /**< Id of molecules.*/
    char *name[2]; /**< Molecule names.*/
    char name_buf[2*CIA_NAME_LEN]; /**< Buffer for molecule names.*/
    fp_t *cross_section; /**< Cross-section [cm^4].*/
    uint64_t num_wpoints; /**< Spectral grid size.*/
    Device_t device; /**< Device id.*/
} CollisionInducedAbsorption_t;


#endif
