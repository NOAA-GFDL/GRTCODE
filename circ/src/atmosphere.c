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
#include "gas_optics.h"
#include "grtcode_utilities.h"
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
#define put_att(a, b, c, d, e, f) nc_catch(nc_put_att_float(a, b, c, d, e ,f))
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#define put_var(a, b, c, d, e) nc_catch(nc_put_vara_double(a, b, c, d, e))
#define put_att(a, b, c, d, e, f) nc_catch(nc_put_att_double(a, b, c, d, e ,f))
#endif
#define reset(start, count) { \
    start[0] = 0; \
    count[0] = 1;}
#define toppmv 1.e6
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/*Read in ppmv values and interpolate from pressure layers to pressure levels.*/
static void get_ppmv(int const ncid, /*Netcdf file id.*/
                     char const * const name, /*Variable name.*/
                     int const level_lower_bound, /*Level lower bound index.*/
                     int const num_levels, /*Number of levels.*/
                     fp_t const * const level_pressure, /*Level pressure (level).*/
                     fp_t const * const layer_pressure, /*Layer pressure (layer).*/
                     fp_t ** const xmol /*Molecule ppmv (level).*/
                    )
{
    int const num_layers = num_levels - 1;
    fp_t *buffer;
    alloc(buffer, num_layers, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, name, &varid));
    size_t start[1];
    size_t count[1];
    reset(start, count);
    start[0] = level_lower_bound;
    count[0] = num_layers;
    get_var(ncid, varid, start, count, buffer);
    fp_t *mol;
    alloc(mol, num_levels, fp_t *);
    mol[0] = buffer[0];
    int j;
    for (j=1; j<num_layers; ++j)
    {
        fp_t const *x = &(layer_pressure[j]);
        fp_t const *y = &(buffer[j]);
        fp_t m = (*y - *(y-1))/(*x - *(x-1));
        fp_t b = *y - (*x)*m;
        mol[j] = m*level_pressure[j] + b;
    }
    mol[num_layers] = buffer[num_layers-1];
    free(buffer);
    *xmol = mol;
    return;
}


/*Reserve memory and read in atmospheric data.*/
void create_atmosphere(Atmosphere_t * const atm, char const * const filepath,
                       int const * const molecules, int const num_molecules,
                       Cfc_t const * const cfc, int const num_cfcs,
                       int const * const cias, int const num_cias)
{
    int ncid;
    nc_catch(nc_open(filepath, NC_NOWRITE, &ncid));

    if (atm->Z < 0)
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "levels", &dimid));
        size_t num_levels;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_levels));
        atm->Z = (int)num_levels - 1;
    }
    atm->num_levels = atm->Z - atm->z + 1;
    atm->num_layers = atm->num_levels - 1;

    alloc(atm->level_pressure, atm->num_levels, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "level_pressure", &varid));
    size_t start[1];
    size_t count[1];
    reset(start, count);
    start[0] = atm->z;
    count[0] = atm->num_levels;
    get_var(ncid, varid, start, count, atm->level_pressure);

    alloc(atm->layer_pressure, atm->num_layers, fp_t *);
    nc_catch(nc_inq_varid(ncid, "layer_pressure", &varid));
    reset(start, count);
    start[0] = atm->z;
    count[0] = atm->num_layers;
    get_var(ncid, varid, start, count, atm->layer_pressure);

    alloc(atm->level_temperature, atm->num_levels, fp_t *);
    nc_catch(nc_inq_varid(ncid, "level_temperature", &varid));
    reset(start, count);
    start[0] = atm->z;
    count[0] = atm->num_levels;
    get_var(ncid, varid, start, count, atm->level_temperature);

    alloc(atm->layer_temperature, atm->num_layers, fp_t *);
    nc_catch(nc_inq_varid(ncid, "layer_temperature", &varid));
    reset(start, count);
    start[0] = atm->z;
    count[0] = atm->num_layers;
    get_var(ncid, varid, start, count, atm->layer_temperature);

    nc_catch(nc_inq_varid(ncid, "surface_temperature", &varid));
    reset(start, count);
    start[0] = 0;
    count[0] = 1;
    get_var(ncid, varid, start, count, &(atm->surface_temperature));

    nc_catch(nc_inq_varid(ncid, "solar_zenith_angle", &varid));
    reset(start, count);
    start[0] = 0;
    count[0] = 1;
    get_var(ncid, varid, start, count, &(atm->solar_zenith_angle));
    atm->solar_zenith_angle = cos(2.*M_PI*atm->solar_zenith_angle/360.);

    nc_catch(nc_inq_varid(ncid, "toa_solar_irradiance", &varid));
    reset(start, count);
    start[0] = 0;
    count[0] = 1;
    get_var(ncid, varid, start, count, &(atm->total_solar_irradiance));
    atm->total_solar_irradiance /= atm->solar_zenith_angle;

    alloc(atm->surface_albedo, atm->grid.n, fp_t *);
    if (atm->alpha >= 0.)
    {
        uint64_t j;
        for (j=0; j<atm->grid.n; ++j)
        {
            atm->surface_albedo[j] = atm->alpha;
        }
    }
    else
    {
        /*Read in the surface albedo and interpolate to the spectral grid.*/
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "wavenumber", &dimid));
        size_t num_wavenumber;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_wavenumber));
        fp_t *w;
        alloc(w, num_wavenumber, fp_t *);
        nc_catch(nc_inq_varid(ncid, "wavenumber", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = num_wavenumber;
        get_var(ncid, varid, start, count, w);
        fp_t *buffer;
        alloc(buffer, num_wavenumber, fp_t *);
        nc_catch(nc_inq_varid(ncid, "surface_albedo", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = num_wavenumber;
        get_var(ncid, varid, start, count, buffer);
        interpolate_to_grid(atm->grid, w, buffer, num_wavenumber, atm->surface_albedo,
                            linear_sample, NULL);
        free(w);
        free(buffer);
    }

    /*CIRC cases assume the surface emissivity is 1.*/
    alloc(atm->surface_emissivity, atm->grid.n, fp_t *);
    uint64_t j;
    for (j=0; j<atm->grid.n; ++j)
    {
        atm->surface_emissivity[j] = 1.;
    }

    char *molecule_names[32];
    molecule_names[CH4] = "CH4_abundance";
    molecule_names[CO] = "CO_abundance";
    molecule_names[CO2] = "CO2_abundance";
    molecule_names[H2O] = "H2O_abundance";
    molecule_names[N2O] = "N2O_abundance";
    molecule_names[O2] = "O2_abundance";
    molecule_names[O3] = "O3_abundance";
    alloc(atm->ppmv, num_molecules, fp_t **);
    atm->num_molecules = num_molecules;
    int i;
    for (i=0; i<num_molecules; ++i)
    {
        fp_t *p;
        get_ppmv(ncid, molecule_names[molecules[i]], atm->z, atm->num_levels,
                 atm->level_pressure, atm->layer_pressure, &p);
        int j;
        for (j=0; j<atm->num_levels; ++j)
        {
            p[j] *= toppmv;
        }
        atm->ppmv[i] = p;
    }

    char *cfc_names[32];
    cfc_names[CFC11] = "CFC11_abundance";
    cfc_names[CFC12] = "CFC12_abundance";
    cfc_names[CCl4] = "CCl4_abundance";
    alloc(atm->cfc_ppmv, num_cfcs, fp_t **);
    atm->num_cfcs = num_cfcs;
    for (i=0; i<num_cfcs; ++i)
    {
        fp_t *p;
        get_ppmv(ncid, cfc_names[cfc[i].id], atm->z, atm->num_levels,
                 atm->level_pressure, atm->layer_pressure, &p);
        int j;
        for (j=0; j<atm->num_levels; ++j)
        {
            p[j] *= toppmv;
        }
        atm->cfc_ppmv[i] = p;
    }

    char *cia_names[32];
    cia_names[CIA_O2] = "O2_abundance";
    alloc(atm->cia_ppmv, num_cias, fp_t **);
    atm->num_cias = num_cias;
    for (i=0; i<num_cias; ++i)
    {
        fp_t *p;
        if (cias[i] == CIA_N2)
        {
            alloc(p, atm->num_levels, fp_t *)
            int j;
            for (j=0; j<atm->num_levels; ++j)
            {
                p[j] = 0.781;
            }
        }
        else
        {
            get_ppmv(ncid, cia_names[cias[i]], atm->z, atm->num_levels,
                     atm->level_pressure, atm->layer_pressure, &p);
        }
        int j;
        for (j=0; j<atm->num_levels; ++j)
        {
            p[j] *= toppmv;
        }
        atm->cia_ppmv[i] = p;
    }

    if (!atm->clean)
    {
        /*Get aerosol optical properties.*/
        fp_t alpha;
        nc_catch(nc_inq_varid(ncid, "angstrom_exponent", &varid));
        reset(start,count);
        start[0] = 0;
        count[0] = 1;
        get_var(ncid, varid, start, count, &alpha);
        fp_t *buffer;
        alloc(buffer, atm->num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "aerosol_optical_depth_at_1_micron", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = atm->num_layers;
        get_var(ncid, varid, start, count, buffer);
        alloc(atm->aerosol_optical_depth, atm->num_layers*atm->grid.n, fp_t *);
        fp_t const cmtomicron = 10000.;
        for (i=0; i<atm->num_layers; ++i)
        {
            for (j=0; j<atm->grid.n; ++j)
            {
                fp_t lambda = cmtomicron/(atm->grid.w0 + j*atm->grid.dw);
                atm->aerosol_optical_depth[i*atm->grid.n+j] = buffer[i]*pow(lambda, -1.*alpha);
            }
        }
        nc_catch(nc_inq_varid(ncid, "aerosol_single_scatter_albedo", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = atm->num_layers;
        get_var(ncid, varid, start, count, buffer);
        alloc(atm->aerosol_single_scatter_albedo, atm->num_layers*atm->grid.n, fp_t *);
        for (i=0; i<atm->num_layers; ++i)
        {
            for (j=0; j<atm->grid.n; ++j)
            {
                atm->aerosol_single_scatter_albedo[i*atm->grid.n+j] = buffer[i];
            }
        }
        nc_catch(nc_inq_varid(ncid, "aerosol_asymmetry_factor", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = atm->num_layers;
        get_var(ncid, varid, start, count, buffer);
        alloc(atm->aerosol_asymmetry_factor, atm->num_layers*atm->grid.n, fp_t *);
        for (i=0; i<atm->num_layers; ++i)
        {
            for (j=0; j<atm->grid.n; ++j)
            {
                atm->aerosol_asymmetry_factor[i*atm->grid.n+j] = buffer[i];
            }
        }
        free(buffer);
    }

    if (!atm->clear)
    {
        /*Get cloud properties.*/
        alloc(atm->liquid_water_path, atm->num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "liquid_water_path", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = atm->num_layers;
        get_var(ncid, varid, start, count, atm->liquid_water_path);
        alloc(atm->liquid_water_droplet_radius, atm->num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "liquid_water_effective_particle_size", &varid));
        reset(start, count);
        start[0] = 0;
        count[0] = atm->num_layers;
        get_var(ncid, varid, start, count, atm->liquid_water_droplet_radius);
    }

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
    free(atm->surface_albedo);
    free(atm->surface_emissivity);
    if (!atm->clean)
    {
        free(atm->aerosol_optical_depth);
        free(atm->aerosol_single_scatter_albedo);
        free(atm->aerosol_asymmetry_factor);
    }
    if (!atm->clear)
    {
        free(atm->liquid_water_path);
        free(atm->liquid_water_droplet_radius);
    }
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
static void add_flux_variable(Output_t * const o, /*Output object.*/
                              char const * const name, /*Variable name.*/
                              char const * const standard_name, /*Variable standard name.*/
                              int const index /*Variable index.*/
                             )
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int varid;
    nc_catch(nc_def_var(o->ncid, name, type, NUM_DIMS, o->dimid, &varid));
    char *unit = "W m-2";
    nc_catch(nc_put_att_text(o->ncid, varid, "units", strlen(unit), unit));
    nc_catch(nc_put_att_text(o->ncid, varid, "standard_name", strlen(standard_name),
                             standard_name));
    o->varid[index] = varid;
    return;
}


/*Create an output file and write metadata.*/
Output_t create_flux_file(char const * const filepath, Atmosphere_t const * const atm)
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
    float zero = 0;
#else
    nc_type type = NC_DOUBLE;
    double zero = 0;
#endif
    Output_t o;
    nc_catch(nc_create(filepath, NC_NETCDF4, &o.ncid));
    nc_catch(nc_def_dim(o.ncid, "level", atm->num_levels, &(o.dimid[LEVEL])));
    add_flux_variable(&o, "rlu", "upwelling_longwave_flux_in_air", RLU);
    add_flux_variable(&o, "rld", "downwelling_longwave_flux_in_air", RLD);
    add_flux_variable(&o, "rsu", "upwelling_shortwave_flux_in_air", RSU);
    put_att(o.ncid, o.varid[RSU], "_FillValue", type, 1, &zero);
    add_flux_variable(&o, "rsd", "downwelling_shortwave_flux_in_air", RSD);
    put_att(o.ncid, o.varid[RSD], "_FillValue", type, 1, &zero);
    return o;
}


/*Close output file.*/
void close_flux_file(Output_t const * const o)
{
    nc_catch(nc_close(o->ncid));
    return;
}


/*Write fluxes to the output file.*/
void write_fluxes(Output_t const * const o, int const index, fp_t const * const flux)
{
    size_t len;
    nc_catch(nc_inq_dimlen(o->ncid, o->dimid[LEVEL], &len));
    size_t start[1];
    size_t count[1];
    reset(start, count);
    start[0] = 0;
    count[0] = len;
    put_var(o->ncid, o->varid[index], start, count, flux);
    return;
}
