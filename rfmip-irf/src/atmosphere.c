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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atmosphere.h"
#include "cfcs.h"
#include "collision_induced_absorption.h"
#include "grtcode_utilities.h"
#include "molecules.h"
#include "netcdf.h"


#define alloc(p, s, t) {p = (t)malloc(sizeof(*p)*s);}
#define nc_catch(e) { \
    int e_ = e; \
    if (e_ != NC_NOERR) {\
        fprintf(stderr, "[%s: %d] %s\n", __FILE__, __LINE__, nc_strerror(e_)); \
        exit(EXIT_FAILURE); \
    }}
#ifdef SINGLE_PRECISION
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_float(a, b, c, d, e))
#define put_var(a, b, c, d, e) nc_catch(nc_put_vara_float(a, b, c, d, e))
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#define put_var(a, b, c, d, e) nc_catch(nc_put_vara_double(a, b, c, d, e))
#endif
#define reset(start, count) { \
    start[0] = 0; start[1] = 0; start[2] = 0; start[3] = 0; \
    count[0] = 1; count[1] = 1; count[2] = 1; count[3] = 1;}
#define Patomb 0.01
#define toppmv 1.e6
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/*Read in ppmv values and interpolate from pressure layers to pressure levels.*/
static void get_ppmv(int const ncid, /*Netcdf file id.*/
                     char const * const name, /*Variable name.*/
                     int const experiment, /*Experiment index.*/
                     int const column_lower_bound, /*Column lower bound index.*/
                     int const num_columns, /*Number of columns.*/
                     int const level_lower_bound, /*Level lower bound index.*/
                     int const num_levels, /*Number of levels.*/
                     fp_t const * const level_pressure, /*Level pressure (column, level).*/
                     fp_t const * const layer_pressure, /*Layer pressure (column, layer).*/
                     fp_t ** const xmol /*Molecule ppmv (column, level).*/
                    )
{
    int const num_layers = num_levels - 1;
    fp_t *buffer;
    alloc(buffer, num_columns*num_layers, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, name, &varid));
    char unit[1024];
    memset(unit, '\0', 1024);
    nc_catch(nc_get_att_text(ncid, varid, "units", unit));
    char *end;
    fp_t unit_val = strtod(unit, &end);
    size_t start[4];
    size_t count[4];
    reset(start, count);
    start[0] = experiment;
    start[1] = column_lower_bound;
    start[2] = level_lower_bound;
    count[1] = num_columns;
    count[2] = num_layers;
    get_var(ncid, varid, start, count, buffer);
    fp_t *mol;
    alloc(mol, num_columns*num_levels, fp_t *);
    int i;
    for (i=0; i<num_columns; ++i)
    {
        int o = i*num_levels;
        int o2 = i*num_layers;
        mol[o] = buffer[o2]*unit_val;
        int j;
        for (j=1; j<num_layers; ++j)
        {
            fp_t const *x = &(layer_pressure[o2+j]);
            fp_t const *y = &(buffer[o2+j]);
            fp_t m = (*y - *(y-1))/(*x - *(x-1));
            fp_t b = *y - (*x)*m;
            mol[o+j] = (m*level_pressure[o+j] + b)*unit_val;
        }
        mol[o+num_layers] = buffer[o2+num_layers-1]*unit_val;
    }
    free(buffer);
    *xmol = mol;
    return;
}


/*Read in "global mean" ppmv values.*/
static void get_gm_ppmv(int const ncid, /*Netcdf file id.*/
                        char const * const name, /*Variable name.*/
                        int const experiment, /*Experiment index.*/
                        int const num_columns, /*Number of columns.*/
                        int const num_levels, /*Number of levels.*/
                        fp_t ** const xmol /*Molecule ppmv (column, level).*/
                       )
{
    int varid;
    nc_catch(nc_inq_varid(ncid, name, &varid));
    char unit[1024];
    memset(unit, '\0', 1024);
    nc_catch(nc_get_att_text(ncid, varid, "units", unit));
    char *end;
    fp_t unit_val = strtod(unit, &end);
    size_t start[4];
    size_t count[4];
    reset(start, count);
    start[0] = experiment;
    fp_t buffer;
    get_var(ncid, varid, start, count, &buffer);
    fp_t *mol;
    alloc(mol, num_columns*num_levels, fp_t *);
    int i;
    for (i=0; i<num_columns*num_levels; ++i)
    {
        mol[i] = buffer*unit_val;
    }
    *xmol = mol;
    return;
}


/*Reserve memory and read in atmospheric data.*/
void create_atmosphere(Atmosphere_t * const atm, char const * const filepath,
                       int const experiment, int const * const molecules,
                       int const num_molecules, Cfc_t const * const cfc,
                       int const num_cfcs, int const * const cias, int const num_cias)
{
    int ncid;
    nc_catch(nc_open(filepath, NC_NOWRITE, &ncid));

    alloc(atm->level_pressure, atm->num_columns*atm->num_levels, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "pres_level", &varid));
    size_t start[4];
    size_t count[4];
    reset(start, count);
    start[0] = atm->x;
    start[1] = atm->z;
    count[0] = atm->num_columns;
    count[1] = atm->num_levels;
    get_var(ncid, varid, start, count, atm->level_pressure);
    int i;
    for (i=0; i<(atm->num_columns*atm->num_levels); ++i)
    {
        atm->level_pressure[i] *= Patomb;
    }

    alloc(atm->layer_pressure, atm->num_columns*atm->num_layers, fp_t *);
    nc_catch(nc_inq_varid(ncid, "pres_layer", &varid));
    reset(start, count);
    start[0] = atm->x;
    start[1] = atm->z;
    count[0] = atm->num_columns;
    count[1] = atm->num_layers;
    get_var(ncid, varid, start, count, atm->layer_pressure);
    for (i=0; i<(atm->num_columns*atm->num_layers); ++i)
    {
        atm->layer_pressure[i] *= Patomb;
    }

    alloc(atm->level_temperature, atm->num_columns*atm->num_levels, fp_t *);
    nc_catch(nc_inq_varid(ncid, "temp_level", &varid));
    reset(start, count);
    start[0] = experiment;
    start[1] = atm->x;
    start[2] = atm->z;
    count[1] = atm->num_columns;
    count[2] = atm->num_levels;
    get_var(ncid, varid, start, count, atm->level_temperature);

    alloc(atm->layer_temperature, atm->num_columns*atm->num_layers, fp_t *);
    nc_catch(nc_inq_varid(ncid, "temp_layer", &varid));
    reset(start, count);
    start[0] = experiment;
    start[1] = atm->x;
    start[2] = atm->z;
    count[1] = atm->num_columns;
    count[2] = atm->num_layers;
    get_var(ncid, varid, start, count, atm->layer_temperature);

    alloc(atm->surface_temperature, atm->num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_temperature", &varid));
    reset(start, count);
    start[0] = experiment;
    start[1] = atm->x;
    count[1] = atm->num_columns;
    get_var(ncid, varid, start, count, atm->surface_temperature);

    alloc(atm->total_solar_irradiance, atm->num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "total_solar_irradiance", &varid));
    reset(start, count);
    start[0] = atm->x;
    count[0] = atm->num_columns;
    get_var(ncid, varid, start, count, atm->total_solar_irradiance);

    alloc(atm->solar_zenith_angle, atm->num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "solar_zenith_angle", &varid));
    reset(start, count);
    start[0] = atm->x;
    count[0] = atm->num_columns;
    get_var(ncid, varid, start, count, atm->solar_zenith_angle);
    for (i=0; i<atm->num_columns; ++i)
    {
        atm->solar_zenith_angle[i] = cos(2.*M_PI*atm->solar_zenith_angle[i]/360.);
    }

    fp_t *buffer;
    alloc(buffer, atm->num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_albedo", &varid));
    reset(start, count);
    start[0] = atm->x;
    count[0] = atm->num_columns;
    get_var(ncid, varid, start, count, buffer);
    alloc(atm->surface_albedo, atm->num_columns*atm->num_sw_wavenumber, fp_t *);
    for (i=0; i<atm->num_columns; ++i)
    {
        uint64_t j;
        for (j=0; j<(atm->num_sw_wavenumber); ++j)
        {
            atm->surface_albedo[i*atm->num_sw_wavenumber+j] = buffer[i];
        }
    }

    nc_catch(nc_inq_varid(ncid, "surface_emissivity", &varid));
    reset(start, count);
    start[0] = atm->x;
    count[0] = atm->num_columns;
    get_var(ncid, varid, start, count, buffer);
    alloc(atm->surface_emissivity, atm->num_columns*atm->num_lw_wavenumber, fp_t *);
    for (i=0; i<atm->num_columns; ++i)
    {
        uint64_t j;
        for (j=0; j<(atm->num_lw_wavenumber); ++j)
        {
            atm->surface_emissivity[i*atm->num_lw_wavenumber+j] = buffer[i];
        }
    }
    free(buffer);

    char *molecule_names[32];
    molecule_names[CH4] = "methane_GM";
    molecule_names[CO] = "carbon_monoxide_GM";
    molecule_names[CO2] = "carbon_dioxide_GM";
    molecule_names[H2O] = "water_vapor";
    molecule_names[N2O] = "nitrous_oxide_GM";
    molecule_names[O2] = "oxygen_GM";
    molecule_names[O3] = "ozone";
    fp_t *xh2o;
    get_ppmv(ncid, molecule_names[H2O], experiment, atm->x, atm->num_columns, atm->z,
             atm->num_levels, atm->level_pressure, atm->layer_pressure, &xh2o);
    alloc(atm->ppmv, num_molecules, fp_t **);
    atm->num_molecules = num_molecules;
    for (i=0; i<num_molecules; ++i)
    {
        fp_t *p;
        if (molecules[i] == H2O)
        {
            alloc(p, atm->num_columns*atm->num_levels, fp_t *);
            memcpy(p, xh2o, sizeof(*p)*atm->num_columns*atm->num_levels);
        }
        else if (molecules[i] == O3)
        {
            get_ppmv(ncid, molecule_names[O3], experiment, atm->x, atm->num_columns, atm->z,
                     atm->num_levels, atm->level_pressure, atm->layer_pressure, &p);
        }
        else
        {
            get_gm_ppmv(ncid, molecule_names[molecules[i]], experiment, atm->num_columns,
                        atm->num_levels, &p);
        }
        int j;
        for (j=0; j<atm->num_columns*atm->num_levels; ++j)
        {
            p[j] *= toppmv/(1. + xh2o[j]);
        }
        atm->ppmv[i] = p;
    }

    alloc(atm->cfc_ppmv, num_cfcs, fp_t **);
    atm->num_cfcs = num_cfcs;
    for (i=0; i<num_cfcs; ++i)
    {
        char *cfc_name;
        switch (cfc[i].id)
        {
            case CFC11:
                if (cfc[i].use_equivalent_ppmv)
                {
                    cfc_name = "cfc11eq_GM";
                }
                else
                {
                    cfc_name = "cfc11_GM";
                }
                break;
            case CFC12:
                if (cfc[i].use_equivalent_ppmv)
                {
                    cfc_name = "cfc12eq_GM";
                }
                else
                {
                    cfc_name = "cfc12_GM";
                }
                break;
            case CFC113:
                cfc_name = "cfc113_GM";
                break;
            case CFC114:
                cfc_name = "cfc114_GM";
                break;
            case CFC115:
                cfc_name = "cfc115_GM";
                break;
            case HCFC22:
                cfc_name = "hcfc22_GM";
                break;
            case HCFC141b:
                cfc_name = "hcfc141b_GM";
                break;
            case HCFC142b:
                cfc_name = "hcfc142b_GM";
                break;
            case HFC23:
                cfc_name = "hfc23_GM";
                break;
            case HFC125:
                cfc_name = "hfc125_GM";
                break;
            case HFC134a:
                if (cfc[i].use_equivalent_ppmv)
                {
                    cfc_name = "hfc134aeq_GM";
                }
                else
                {
                    cfc_name = "hfc134a_GM";
                }
                break;
            case HFC143a:
                cfc_name = "hfc143a_GM";
                break;
            case HFC152a:
                cfc_name = "hfc152a_GM";
                break;
            case HFC227ea:
                cfc_name = "hfc227ea_GM";
                break;
            case HFC245fa:
                cfc_name = "hfc245fa_GM";
                break;
            case CCl4:
                cfc_name = "carbon_tetrachloride_GM";
                break;
            case C2F6:
                cfc_name = "c2f6_GM";
                break;
            case CF4:
                cfc_name = "cf4_GM";
                break;
            case CH2Cl2:
                cfc_name = "ch2cl2_GM";
                break;
            case NF3:
                cfc_name = "nf3_GM";
                break;
            case SF6:
                cfc_name = "sf6_GM";
                break;
            default:
                fprintf(stderr, "[%s: %d] unknown CFC id.\n", __FILE__, __LINE__);
                exit(EXIT_FAILURE);
        }
        fp_t *p;
        get_gm_ppmv(ncid, cfc_name, experiment, atm->num_columns, atm->num_levels, &p);
        int j;
        for (j=0; j<atm->num_columns*atm->num_levels; ++j)
        {
            p[j] *= toppmv/(1. + xh2o[j]);
        }
        atm->cfc_ppmv[i] = p;
    }

    char *cia_names[32];
    cia_names[CIA_N2] = "nitrogen_GM";
    cia_names[CIA_O2] = "oxygen_GM";
    alloc(atm->cia_ppmv, num_cias, fp_t **);
    atm->num_cias = num_cias;
    for (i=0; i<num_cias; ++i)
    {
        fp_t *p;
        get_gm_ppmv(ncid, cia_names[cias[i]], experiment, atm->num_columns,
                    atm->num_levels, &p);
        int j;
        for (j=0; j<atm->num_columns*atm->num_levels; ++j)
        {
            p[j] *= toppmv/(1. + xh2o[j]);
        }
        atm->cia_ppmv[i] = p;
    }

    free(xh2o);
    nc_catch(nc_close(ncid));
    return;
}


/*Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm)
{
    free(atm->level_pressure);
    free(atm->layer_pressure);
    free(atm->level_temperature);
    free(atm->layer_temperature);
    free(atm->surface_temperature);
    free(atm->total_solar_irradiance);
    free(atm->solar_zenith_angle);
    free(atm->surface_albedo);
    free(atm->surface_emissivity);
    int i;
    for (i=0; i<atm->num_molecules; ++i)
    {
        free(atm->ppmv[i]);
    }
    free(atm->ppmv);
    for (i=0; i<atm->num_cfcs; ++i)
    {
        free(atm->cfc_ppmv[i]);
    }
    free(atm->cfc_ppmv);
    for (i=0; i<atm->num_cias; ++i)
    {
        free(atm->cia_ppmv[i]);
    }
    free(atm->cia_ppmv);
    return;
}


/*Add a variable to the output file.*/
static int add_variable(int const ncid, /*File id.*/
                        char const * const name, /*Name.*/
                        char const * const standard_name, /*Standard name.*/
                        char const * const units, /*Units.*/
                        int const * const dimid,
                        int const num_dims
                       )
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int varid;
    nc_catch(nc_def_var(ncid, name, type, num_dims, dimid, &varid));
    nc_catch(nc_put_att_text(ncid, varid, "units", strlen(units), units));
    nc_catch(nc_put_att_text(ncid, varid, "standard_name", strlen(standard_name),
                             standard_name));
    return varid;
}


/*Create an output file and write metadata.*/
Output_t create_flux_file(char const * const filepath, Atmosphere_t const * const atm)
{
    Output_t o;
    nc_catch(nc_create(filepath, NC_NETCDF4, &o.ncid));
    nc_catch(nc_def_dim(o.ncid, "expt", 1, &(o.dimid[EXPT])));
    nc_catch(nc_def_dim(o.ncid, "site", atm->num_columns, &(o.dimid[SITE])));
    nc_catch(nc_def_dim(o.ncid, "level", atm->num_levels, &(o.dimid[LEVEL])));
    int dimid[3] = {o.dimid[EXPT], o.dimid[SITE], o.dimid[LEVEL]};
    o.varid[RLU] = add_variable(o.ncid, "rlu", "upwelling_longwave_flux_in_air",
                                "W m-2", dimid, 3);
    o.varid[RLD] = add_variable(o.ncid, "rld", "downwelling_longwave_flux_in_air",
                                "W m-2", dimid, 3);
    o.varid[RSU] = add_variable(o.ncid, "rsu", "upwelling_shortwave_flux_in_air",
                                "W m-2", dimid, 3);
    o.varid[RSD] = add_variable(o.ncid, "rsd", "downwelling_shortwave_flux_in_air",
                                "W m-2", dimid, 3);
    return o;
}


/*Close output file.*/
void close_flux_file(Output_t const * const o)
{
    nc_catch(nc_close(o->ncid));
    return;
}


/*Write a column of fluxes to the output file.*/
void write_output(Output_t const * const o, int const index, fp_t const * const data,
                  size_t const * const start, size_t const * const count)
{
    put_var(o->ncid, o->varid[index], start, count, data);
    return;
}
