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

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "debug.h"
#include "grtcode_utilities.h"
#include "molecules.h"
#include "molecules-internal.h"
#include "parse_HITRAN_file-internal.h"


/*Initialize a molecule.  Read in its line parameters from the input HITRAN database file.*/
int create_molecule(Molecule_t * const mol, int const id, char const * const hitran_path,
                    double const w0, double const wn, int const num_layers,
                    Device_t const device)
{
    not_null(mol);
    mol->id = id;
    switch(mol->id)
    {
        case H2O:
            snprintf(mol->name, MOL_NAME_LEN, "H2O");
            mol->mass = 18.010565f;
            mol->num_isotopologues = 9;
            break;
        case CO2:
            snprintf(mol->name, MOL_NAME_LEN, "CO2");
            mol->mass = 43.98983f;
            mol->num_isotopologues = 13;
            break;
        case O3:
            snprintf(mol->name, MOL_NAME_LEN, "O3");
            mol->mass = 47.984745f;
            mol->num_isotopologues = 18;
            break;
        case N2O:
            snprintf(mol->name, MOL_NAME_LEN, "N2O");
            mol->mass = 44.001062f;
            mol->num_isotopologues = 5;
            break;
        case CO:
            snprintf(mol->name, MOL_NAME_LEN, "CO");
            mol->mass = 27.994915f;
            mol->num_isotopologues = 9;
            break;
        case CH4:
            snprintf(mol->name, MOL_NAME_LEN, "CH4");
            mol->mass = 16.0313f;
            mol->num_isotopologues = 4;
            break;
        case O2:
            snprintf(mol->name, MOL_NAME_LEN, "O2");
            mol->mass = 31.98983f;
            mol->num_isotopologues = 6;
            break;
        case NO:
            snprintf(mol->name, MOL_NAME_LEN, "NO");
            mol->mass = 29.997989f;
            mol->num_isotopologues = 3;
            break;
        case SO2:
            snprintf(mol->name, MOL_NAME_LEN, "SO2");
            mol->mass = 63.961901f;
            mol->num_isotopologues = 2;
            break;
        case NO2:
            snprintf(mol->name, MOL_NAME_LEN, "NO2");
            mol->mass = 45.992904f;
            mol->num_isotopologues = 1;
            break;
        case NH3:
            snprintf(mol->name, MOL_NAME_LEN, "NH3");
            mol->mass = 17.026549f;
            mol->num_isotopologues = 2;
            break;
        case HNO3:
            snprintf(mol->name, MOL_NAME_LEN, "HNO3");
            mol->mass = 62.995644f;
            mol->num_isotopologues = 2;
            break;
        case OH:
            snprintf(mol->name, MOL_NAME_LEN, "OH");
            mol->mass = 17.00274f;
            mol->num_isotopologues = 3;
            break;
        case HF:
            snprintf(mol->name, MOL_NAME_LEN, "HF");
            mol->mass = 20.006229f;
            mol->num_isotopologues = 2;
            break;
        case HCl:
            snprintf(mol->name, MOL_NAME_LEN, "HCl");
            mol->mass = 35.976678f;
            mol->num_isotopologues = 4;
            break;
        case HBr:
            snprintf(mol->name, MOL_NAME_LEN, "HBr");
            mol->mass = 79.92616f;
            mol->num_isotopologues = 4;
            break;
        case HI:
            snprintf(mol->name, MOL_NAME_LEN, "HI");
            mol->mass = 127.912297f;
            mol->num_isotopologues = 2;
            break;
        case ClO:
            snprintf(mol->name, MOL_NAME_LEN, "ClO");
            mol->mass = 50.963768f;
            mol->num_isotopologues = 2;
            break;
        case OCS:
            snprintf(mol->name, MOL_NAME_LEN, "OCS");
            mol->mass = 59.966986f;
            mol->num_isotopologues = 5;
            break;
        case H2CO:
            snprintf(mol->name, MOL_NAME_LEN, "H2CO");
            mol->mass = 30.010565f;
            mol->num_isotopologues = 3;
            break;
        case HOCl:
            snprintf(mol->name, MOL_NAME_LEN, "HOCl");
            mol->mass = 51.971593f;
            mol->num_isotopologues = 2;
            break;
        case N2:
            snprintf(mol->name, MOL_NAME_LEN, "N2");
            mol->mass = 28.006148f;
            mol->num_isotopologues = 3;
            break;
        case HCN:
            snprintf(mol->name, MOL_NAME_LEN, "HCN");
            mol->mass = 27.010899f;
            mol->num_isotopologues = 3;
            break;
        case CH3Cl:
            snprintf(mol->name, MOL_NAME_LEN, "CH3Cl");
            mol->mass = 49.992328f;
            mol->num_isotopologues = 2;
            break;
        case H2O2:
            snprintf(mol->name, MOL_NAME_LEN, "H2O2");
            mol->mass = 34.00548f;
            mol->num_isotopologues = 1;
            break;
        case C2H2:
            snprintf(mol->name, MOL_NAME_LEN, "C2H2");
            mol->mass = 26.01565f;
            mol->num_isotopologues = 3;
            break;
        case C2H6:
            snprintf(mol->name, MOL_NAME_LEN, "C2H6");
            mol->mass = 30.04695f;
            mol->num_isotopologues = 2;
            break;
        case PH3:
            snprintf(mol->name, MOL_NAME_LEN, "PH3");
            mol->mass = 33.997238f;
            mol->num_isotopologues = 1;
            break;
        case COF2:
            snprintf(mol->name, MOL_NAME_LEN, "COF2");
            mol->mass = 65.991722f;
            mol->num_isotopologues = 2;
            break;
        case SF6_MOL:
            snprintf(mol->name, MOL_NAME_LEN, "SF6");
            mol->mass = 145.962492f;
            mol->num_isotopologues = 1;
            break;
        case H2S:
            snprintf(mol->name, MOL_NAME_LEN, "H2S");
            mol->mass = 33.987721f;
            mol->num_isotopologues = 3;
            break;
        case HCOOH:
            snprintf(mol->name, MOL_NAME_LEN, "HCOOH");
            mol->mass = 46.00548f;
            mol->num_isotopologues = 1;
            break;
        case HO2:
            snprintf(mol->name, MOL_NAME_LEN, "HO2");
            mol->mass = 32.997655f;
            mol->num_isotopologues = 1;
            break;
        case O:
            snprintf(mol->name, MOL_NAME_LEN, "O");
            mol->mass = 15.994915f;
            mol->num_isotopologues = 0;
            break;
        case ClONO2:
            snprintf(mol->name, MOL_NAME_LEN, "ClONO2");
            mol->mass = 96.956672f;
            mol->num_isotopologues = 2;
            break;
        case NOp:
            snprintf(mol->name, MOL_NAME_LEN, "NO+");
            mol->mass = 29.997989f;
            mol->num_isotopologues = 1;
            break;
        case HOBr:
            snprintf(mol->name, MOL_NAME_LEN, "HOBr");
            mol->mass = 95.921076f;
            mol->num_isotopologues = 2;
            break;
        case C2H4:
            snprintf(mol->name, MOL_NAME_LEN, "C2H4");
            mol->mass = 28.0313f;
            mol->num_isotopologues = 2;
            break;
        case CH3OH:
            snprintf(mol->name, MOL_NAME_LEN, "CH3OH");
            mol->mass = 32.026215f;
            mol->num_isotopologues = 1;
            break;
        case CH3Br:
            snprintf(mol->name, MOL_NAME_LEN, "CH3Br");
            mol->mass = 93.941811f;
            mol->num_isotopologues = 2;
            break;
        case CH3CN:
            snprintf(mol->name, MOL_NAME_LEN, "CH3CN");
            mol->mass = 41.026549f;
            mol->num_isotopologues = 4;
            break;
        case CF4_MOL:
            snprintf(mol->name, MOL_NAME_LEN, "CF4");
            mol->mass = 87.993616f;
            mol->num_isotopologues = 1;
            break;
        case C4H2:
            snprintf(mol->name, MOL_NAME_LEN, "C4H2");
            mol->mass = 50.01565f;
            mol->num_isotopologues = 1;
            break;
        case HC3N:
            snprintf(mol->name, MOL_NAME_LEN, "HC3N");
            mol->mass = 51.010899f;
            mol->num_isotopologues = 6;
            break;
        case H2:
            snprintf(mol->name, MOL_NAME_LEN, "H2");
            mol->mass = 2.01565f;
            mol->num_isotopologues = 2;
            break;
        case CS:
            snprintf(mol->name, MOL_NAME_LEN, "CS");
            mol->mass = 43.971036f;
            mol->num_isotopologues = 4;
            break;
        case SO3:
            snprintf(mol->name, MOL_NAME_LEN, "SO3");
            mol->mass = 79.95682f;
            mol->num_isotopologues = 1;
            break;
        case C2N2:
            snprintf(mol->name, MOL_NAME_LEN, "C2N2");
            mol->mass = 52.006148f;
            mol->num_isotopologues = 2;
            break;
        case COCl2:
            snprintf(mol->name, MOL_NAME_LEN, "COCl2");
            mol->mass = 97.9326199796f;
            mol->num_isotopologues = 2;
            break;
        case SO:
            snprintf(mol->name, MOL_NAME_LEN, "SO");
            mol->mass = 48.0644f;
            mol->num_isotopologues = 3;
            break;
        case C3H4:
            snprintf(mol->name, MOL_NAME_LEN, "C3H4");
            mol->mass = 40.0639f;
            mol->num_isotopologues = 1;
            break;
        case CH3:
            snprintf(mol->name, MOL_NAME_LEN, "CH3");
            mol->mass = 15.035f;
            mol->num_isotopologues = 1;
            break;
        case CS2:
            snprintf(mol->name, MOL_NAME_LEN, "CS2");
            mol->mass = 76.139f;
            mol->num_isotopologues = 4;
            break;
        default:
            {char const *mesg = "unrecognized molecule id %d.";
            raise(GRTCODE_VALUE_ERR, mesg, mol->id);}
    }
    mol->mass /= 6.023e23;
    mol->device = device;
    catch(parse_hitran_file(&(mol->line_params), hitran_path, mol->id,
                            w0, wn, device));
    gmalloc(mol->q, mol->num_isotopologues*num_layers, device);
    return GRTCODE_SUCCESS;
}


/*Free memory stored by the molecule object.*/
int free_molecule(Molecule_t * const mol)
{
    not_null(mol);
    catch(free_line_params(&(mol->line_params)));
    gfree(mol->q, mol->device);
    return GRTCODE_SUCCESS;
}


/*Convert from a molecule id to an array index.*/
int molecule_hash(int const mol_id, int * const hash)
{
    not_null(hash);
    if (mol_id < H2O || mol_id > NUM_MOLS)
    {
        char const *mesg = "unrecognized molecule id %d.";
        raise(GRTCODE_VALUE_ERR, mesg, mol_id);
    }
    *hash = mol_id - 1;
    return GRTCODE_SUCCESS;
}
