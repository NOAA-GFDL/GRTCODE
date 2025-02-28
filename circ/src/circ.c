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
#include "circ.h"
#include "driver.h"
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


enum dimid
{
    LEVEL = 0,
    NUM_DIMS
};


struct Output
{
    int *dimid;
    int ncid;
    int *varid;
};


/*Reserve memory and read in atmospheric data.*/
Atmosphere_t create_atmosphere(Parser_t *parser)
{
    /*Add/parse command line arguments.*/
    snprintf(parser->description, desclen,
             "Calculates the radiative fluxes for the NASA CIRC test cases.");
    add_argument(parser, "input_file", NULL, "Input data file.", NULL);
    int one = 1;
    add_argument(parser, "-a", "--surface-albedo", "Spectrally constant surface albedo.", &one);
    add_argument(parser, "-CCl4", NULL, "CSV file with CCl4 cross sections.", &one);
    add_argument(parser, "-CFC-11", NULL, "CSV file with CFC-11 cross sections.", &one);
    add_argument(parser, "-CFC-12", NULL, "CSV file with CFC-12 cross sections.", &one);
    add_argument(parser, "-CH4", NULL, "Include CH4.", NULL);
    add_argument(parser, "-clean", NULL, "Run without aerosols.", NULL);
    add_argument(parser, "-clear", NULL, "Run without clouds.", NULL);
    add_argument(parser, "-CO", NULL, "Include CO.", NULL);
    add_argument(parser, "-CO2", NULL, "Include CO2.", NULL);
    add_argument(parser, "-H2O", NULL, "Include H2O.", NULL);
    add_argument(parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(parser, "-N2-N2", NULL, "CSV file with N2-N2 collison cross sections", &one);
    add_argument(parser, "-N2O", NULL, "Include N2O.", NULL);
    add_argument(parser, "-O2", NULL, "Include O2.", NULL);
    add_argument(parser, "-O2-N2", NULL, "CSV file with O2-N2 collison cross sections", &one);
    add_argument(parser, "-O2-O2", NULL, "CSV file with O2-O2 collison cross sections", &one);
    add_argument(parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(parser, "-p", "--cloud", "Cloud parameterization.", &one);
    add_argument(parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(*parser);

    /*Open the input file.*/
    char buffer[valuelen];
    get_argument(*parser, "input_file", buffer);
    int ncid;
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));

    /*Determine the number of levels.*/
    Atmosphere_t atm;
    int z = get_argument(*parser, "-z", buffer) ? atoi(buffer) : 0;
    int Z;
    if (get_argument(*parser, "-Z", buffer))
    {
        Z = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "levels", &dimid));
        size_t num_levels;
        nc_catch(nc_inq_dimlen(ncid, dimid, &num_levels));
        Z = (int)num_levels - 1;
    }
    atm.num_levels = Z - z + 1;
    atm.num_layers = atm.num_levels - 1;
    atm.num_columns = 1;
    atm.num_times = 1;

    /*Pressure.*/
    alloc(atm.level_pressure, atm.num_levels*atm.num_columns, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "level_pressure", &varid));
    size_t start = z;
    size_t count = atm.num_levels*atm.num_columns;
    get_var(ncid, varid, &start, &count, atm.level_pressure);
    alloc(atm.layer_pressure, atm.num_layers*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "layer_pressure", &varid));
    count = atm.num_layers*atm.num_columns;
    get_var(ncid, varid, &start, &count, atm.layer_pressure);

    /*Temperature.*/
    alloc(atm.level_temperature, atm.num_levels*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "level_temperature", &varid));
    count = atm.num_levels*atm.num_columns;
    get_var(ncid, varid, &start, &count, atm.level_temperature);
    alloc(atm.layer_temperature, atm.num_layers*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "layer_temperature", &varid));
    count = atm.num_layers*atm.num_columns;
    get_var(ncid, varid, &start, &count, atm.layer_temperature);
    alloc(atm.surface_temperature, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_temperature", &varid));
    start = 0;
    count = atm.num_columns;
    get_var(ncid, varid, &start, &count, atm.surface_temperature);

    /*Solar zenith angle.*/
    alloc(atm.solar_zenith_angle, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "solar_zenith_angle", &varid));
    get_var(ncid, varid, &start, &count, atm.solar_zenith_angle);
    int i;
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.solar_zenith_angle[i] = cos(2.*M_PI*atm.solar_zenith_angle[i]/360.);
    }

    /*Solar irradiance.*/
    alloc(atm.total_solar_irradiance, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "toa_solar_irradiance", &varid));
    get_var(ncid, varid, &start, &count, atm.total_solar_irradiance);
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.total_solar_irradiance[i] /= atm.solar_zenith_angle[i];
    }

    /*Surface albedo.*/
    if (get_argument(*parser, "-a", buffer))
    {
        atm.albedo_grid_size = 2;
        alloc(atm.albedo_grid, atm.albedo_grid_size, fp_t *);
        atm.albedo_grid[0] = -1.;
        atm.albedo_grid[1] = 0.;
        alloc(atm.surface_albedo, atm.albedo_grid_size, fp_t *);
        atm.surface_albedo[0] = (fp_t)(atof(buffer));
        atm.surface_albedo[1] = atm.surface_albedo[0];
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "wavenumber", &dimid));
        nc_catch(nc_inq_dimlen(ncid, dimid, &atm.albedo_grid_size));
        alloc(atm.albedo_grid, atm.albedo_grid_size, fp_t *);
        nc_catch(nc_inq_varid(ncid, "wavenumber", &varid));
        count = atm.albedo_grid_size;
        get_var(ncid, varid, &start, &count, atm.albedo_grid);
        alloc(atm.surface_albedo, atm.albedo_grid_size, fp_t *);
        nc_catch(nc_inq_varid(ncid, "surface_albedo", &varid));
        get_var(ncid, varid, &start, &count, atm.surface_albedo);
    }

    /*Surface emissivity. CIRC cases assume this is 1.*/
    atm.emissivity_grid_size = 2;
    alloc(atm.emissivity_grid, atm.emissivity_grid_size, fp_t *);
    atm.emissivity_grid[0] = -1.;
    atm.emissivity_grid[1] = 0.;
    alloc(atm.surface_emissivity, atm.emissivity_grid_size, fp_t *);
    atm.surface_emissivity[0] = 1.;
    atm.surface_emissivity[1] = atm.surface_emissivity[0];

    /*Molecular abundances.*/
    fp_t const to_ppmv = 1.e6;
    struct MoleculeMeta
    {
        int id;
        char *flag;
        char *name;
    };
    int const num_molecules = 7;
    struct MoleculeMeta molecules[num_molecules] = {{CH4, "-CH4", "CH4_abundance"},
        {CO, "-CO", "CO_abundance"}, {CO2, "-CO2", "CO2_abundance"},
        {H2O, "-H2O", "H2O_abundance"}, {N2O, "-N2O", "N2O_abundance"},
        {O2, "-O2", "O2_abundance"}, {O3, "-O3", "O3_abundance"}};
    alloc(atm.molecules, num_molecules, int *);
    atm.num_molecules = 0;
    alloc(atm.ppmv, num_molecules, fp_t **);
    fp_t *abundance;
    alloc(abundance, atm.num_layers, fp_t *);
    for (i=0; i<num_molecules; ++i)
    {
        if (get_argument(*parser, molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = molecules[i].id;
            alloc(atm.ppmv[atm.num_molecules], atm.num_levels*atm.num_columns, fp_t *);
            nc_catch(nc_inq_varid(ncid, molecules[i].name, &varid));
            start = z;
            count = atm.num_layers;
            get_var(ncid, varid, &start, &count, abundance);
            fp_t *ppmv = atm.ppmv[atm.num_molecules];
            ppmv[0] = abundance[0]*to_ppmv;
            ppmv[atm.num_layers] = abundance[atm.num_layers-1]*to_ppmv;
            int j;
            for (j=1; j<atm.num_layers; ++j)
            {
                ppmv[j] = to_ppmv*(abundance[j] + (abundance[j+1] - abundance[j])*
                          (atm.level_pressure[j] - atm.layer_pressure[j])/
                          (atm.layer_pressure[j+1] - atm.layer_pressure[j]));
            }
            atm.num_molecules++;
        }
    }

    /*Molecular continua.*/
    if (!get_argument(*parser, "-h2o-ctm", atm.h2o_ctm))
    {
        snprintf(atm.h2o_ctm, valuelen, "%s", "none");
    }
    if (!get_argument(*parser, "-o3-ctm", atm.o3_ctm))
    {
        snprintf(atm.o3_ctm, valuelen, "%s", "none");
    }

    /*CFC abundances.*/
    int const num_cfcs = 3;
    struct MoleculeMeta cfcs[num_cfcs] = {{CFC11, "-CFC-11", "CFC11_abundance"},
        {CFC12, "-CFC-12", "CFC12_abundance"}, {CCl4, "-CCl4", "CCl4_abundance"}};
    alloc(atm.cfc, num_cfcs, Cfc_t *);
    atm.num_cfcs = 0;
    alloc(atm.cfc_ppmv, num_cfcs, fp_t **);
    for (i=0; i<num_cfcs; ++i)
    {
        if (get_argument(*parser, cfcs[i].flag, atm.cfc[atm.num_cfcs].path))
        {
            atm.cfc[atm.num_cfcs].id = cfcs[i].id;
            alloc(atm.cfc_ppmv[atm.num_cfcs], atm.num_levels*atm.num_columns, fp_t *);
            nc_catch(nc_inq_varid(ncid, cfcs[i].name, &varid));
            start = z;
            count = atm.num_layers;
            get_var(ncid, varid, &start, &count, abundance);
            fp_t *ppmv = atm.cfc_ppmv[atm.num_cfcs];
            ppmv[0] = abundance[0]*to_ppmv;
            ppmv[atm.num_layers] = abundance[atm.num_layers-1]*to_ppmv;
            int j;
            for (j=1; j<atm.num_layers; ++j)
            {
                ppmv[j] = to_ppmv*(abundance[j] + (abundance[j+1] - abundance[j])*
                          (atm.level_pressure[j] - atm.layer_pressure[j])/
                          (atm.layer_pressure[j+1] - atm.layer_pressure[j]));
            }
            atm.num_cfcs++;
        }
    }

    /*Collision-induced absorption (CIA) abundances.*/
    int const num_cias = 3;
    int const num_cia_species = 2;
    struct CiaMeta
    {
        int species1;
        int species2;
        char *flag;
    };
    struct CiaMeta cias[num_cias] = {{CIA_N2, CIA_N2, "-N2-N2"},
        {CIA_O2, CIA_N2, "-O2-N2"}, {CIA_O2, CIA_O2, "-O2-O2"}};
    alloc(atm.cia, num_cias, Cia_t *);
    atm.num_cias = 0;
    alloc(atm.cia_species, num_cia_species, int *);
    atm.num_cia_species = 0;
    alloc(atm.cia_ppmv, num_cia_species, fp_t **);
    for (i=0; i<num_cias; ++i)
    {
        if (get_argument(*parser, cias[i].flag, atm.cia[atm.num_cias].path))
        {
            atm.cia[atm.num_cias].id[0] = cias[i].species1;
            atm.cia[atm.num_cias].id[1] = cias[i].species2;
            int j;
            for (j=0; j<2; ++j)
            {
                int k;
                for (k=0; k<atm.num_cia_species; ++k)
                {
                    if (atm.cia[atm.num_cias].id[j] == atm.cia_species[k])
                    {
                        break;
                    }
                }
                if (k >= atm.num_cia_species)
                {
                    atm.cia_species[atm.num_cia_species] = atm.cia[atm.num_cias].id[j];
                    alloc(atm.cia_ppmv[atm.num_cia_species], atm.num_levels*atm.num_columns, fp_t *);
                    fp_t *ppmv = atm.cia_ppmv[atm.num_cia_species];
                    if (atm.cia_species[atm.num_cia_species] == CIA_N2)
                    {
                        for (k=0; k<atm.num_levels; ++k)
                        {
                            ppmv[k] = 0.781*to_ppmv;
                        }
                    }
                    else if (atm.cia_species[atm.num_cia_species] == CIA_O2)
                    {
                        nc_catch(nc_inq_varid(ncid, "O2_abundance", &varid));
                        start = z;
                        count = atm.num_layers;
                        get_var(ncid, varid, &start, &count, abundance);
                        ppmv[0] = abundance[0]*to_ppmv;
                        ppmv[atm.num_layers] = abundance[atm.num_layers-1]*to_ppmv;
                        for (k=1; k<atm.num_layers; ++k)
                        {
                            ppmv[k] = to_ppmv*(abundance[k] + (abundance[k+1] - abundance[k])*
                                      (atm.level_pressure[k] - atm.layer_pressure[k])/
                                      (atm.layer_pressure[k+1] - atm.layer_pressure[k]));
                        }
                    }
                    atm.num_cia_species++;
                }
            }
            atm.num_cias++;
        }
    }
    free(abundance);

    /*Aerosols.*/
    atm.clean = get_argument(*parser, "-clean", NULL);
    if (!atm.clean)
    {
        nc_catch(nc_inq_varid(ncid, "angstrom_exponent", &varid));
        fp_t alpha;
        start = 0;
        count = 1;
        get_var(ncid, varid, &start, &count, &alpha);
        atm.aerosol_grid_size = 50000;
        alloc(atm.aerosol_grid, atm.aerosol_grid_size, fp_t *);
        for (i=0; i<atm.aerosol_grid_size; ++i)
        {
            atm.aerosol_grid[i] = (fp_t)(i + 1);
        }
        fp_t *optics;
        alloc(optics, atm.num_layers, fp_t *);
        nc_catch(nc_inq_varid(ncid, "aerosol_optical_depth_at_1_micron", &varid));
        start = z;
        count = atm.num_layers;
        get_var(ncid, varid, &start, &count, optics);
        alloc(atm.aerosol_optical_depth, atm.num_layers*atm.aerosol_grid_size, fp_t *);
        fp_t const cmtomicron = 10000.;
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                fp_t lambda = cmtomicron/(atm.aerosol_grid[j]);
                atm.aerosol_optical_depth[i*atm.aerosol_grid_size+j] = optics[i]*pow(lambda, -1.*alpha);
            }
        }
        nc_catch(nc_inq_varid(ncid, "aerosol_single_scatter_albedo", &varid));
        get_var(ncid, varid, &start, &count, optics);
        alloc(atm.aerosol_single_scatter_albedo, atm.num_layers*atm.aerosol_grid_size, fp_t *);
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                atm.aerosol_single_scatter_albedo[i*atm.aerosol_grid_size+j] = optics[i];
            }
        }
        nc_catch(nc_inq_varid(ncid, "aerosol_asymmetry_factor", &varid));
        get_var(ncid, varid, &start, &count, optics);
        alloc(atm.aerosol_asymmetry_factor, atm.num_layers*atm.aerosol_grid_size, fp_t *);
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                atm.aerosol_asymmetry_factor[i*atm.aerosol_grid_size+j] = optics[i];
            }
        }
        free(optics);
    }

    /*Clouds.*/
    atm.clear = get_argument(*parser, "-clear", NULL);
    if (!atm.clear)
    {
        alloc(atm.liquid_water_content, atm.num_layers*atm.num_columns, fp_t *);
        nc_catch(nc_inq_varid(ncid, "liquid_water_path", &varid));
        start = z;
        count = atm.num_layers;
        get_var(ncid, varid, &start, &count, atm.liquid_water_content);
        alloc(atm.liquid_water_droplet_radius, atm.num_layers*atm.num_columns, fp_t *);
        nc_catch(nc_inq_varid(ncid, "liquid_water_effective_particle_size", &varid));
        start = z;
        count = atm.num_layers;
        get_var(ncid, varid, &start, &count, atm.liquid_water_droplet_radius);
    }

    /*Close input file.*/
    nc_catch(nc_close(ncid));
    return atm;
}


/*Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm)
{
    free(atm->layer_pressure);
    free(atm->level_pressure);
    free(atm->layer_temperature);
    free(atm->level_temperature);
    free(atm->solar_zenith_angle);
    free(atm->surface_temperature);
    free(atm->total_solar_irradiance);
    free(atm->albedo_grid);
    free(atm->surface_albedo);
    free(atm->emissivity_grid);
    free(atm->surface_emissivity);
    if (!atm->clean)
    {
        free(atm->aerosol_grid);
        free(atm->aerosol_optical_depth);
        free(atm->aerosol_single_scatter_albedo);
        free(atm->aerosol_asymmetry_factor);
    }
    if (!atm->clear)
    {
        free(atm->liquid_water_droplet_radius);
        free(atm->liquid_water_content);
    }
    int i;
    for (i=0; i<atm->num_molecules; ++i)
    {
        free(atm->ppmv[i]);
    }
    free(atm->ppmv);
    free(atm->molecules);
    for (i=0; i<atm->num_cfcs; ++i)
    {
        free(atm->cfc_ppmv[i]);
    }
    free(atm->cfc_ppmv);
    free(atm->cfc);
    for (i=0; i<atm->num_cia_species; ++i)
    {
        free(atm->cia_ppmv[i]);
    }
    free(atm->cia_ppmv);
    free(atm->cia);
    free(atm->cia_species);
    return;
}


/*Add a variable to the output file.*/
static void add_flux_variable(Output_t * const o, /*Output object.*/
                              VarId_t const index, /*Variable index.*/
                              char const * const name, /*Variable name.*/
                              char const * const standard_name, /*Variable standard name.*/
                              fp_t const * const fill_value /*Fill value.*/
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
    if (fill_value != NULL)
    {
#ifdef SINGLE_PRESCISION
        nc_catch(nc_put_att_float(o->ncid, varid, "_FillValue", type, 1, fill_value));
#else
        nc_catch(nc_put_att_double(o->ncid, varid, "_FillValue", type, 1, fill_value));
#endif
    }
    o->varid[index] = varid;
    return;
}


/*Create an output file and write metadata.*/
void create_flux_file(Output_t **output, char const * const filepath,
                      Atmosphere_t const * const atm, SpectralGrid_t const * const lw_grid,
                      SpectralGrid_t const * const sw_grid)
{
    Output_t *file = (Output_t *)malloc(sizeof(*file));
    file->dimid = (int *)malloc(sizeof(*(file->dimid))*NUM_DIMS);
    file->varid = (int *)malloc(sizeof(*(file->varid))*NUM_VARS);
    nc_catch(nc_create(filepath, NC_NETCDF4, &(file->ncid)));
    nc_catch(nc_def_dim(file->ncid, "level", atm->num_levels, &(file->dimid[LEVEL])));
    add_flux_variable(file, RLUAF, "rlu", "upwelling_longwave_flux_in_air", NULL);
    add_flux_variable(file, RLDAF, "rld", "downwelling_longwave_flux_in_air", NULL);
    fp_t const zero = 0;
    add_flux_variable(file, RSU, "rsu", "upwelling_shortwave_flux_in_air", &zero);
    add_flux_variable(file, RSD, "rsd", "downwelling_shortwave_flux_in_air", &zero);
    *output = file;
    return;
}


/*Close output file.*/
void close_flux_file(Output_t * const o)
{
    nc_catch(nc_close(o->ncid));
    free(o->varid);
    free(o->dimid);
    return;
}


/*Write fluxes to the output file.*/
void write_output(Output_t *output, VarId_t id, fp_t *data, int time, int column)
{
    size_t num_levels;
    nc_catch(nc_inq_dimlen(output->ncid, output->dimid[LEVEL], &num_levels));
    size_t start = 0;
    size_t count = num_levels;
#ifdef SINGLE_PRECISION
    nc_catch(nc_put_vara_float(output->ncid, output->varid[id], &start, &count, data))
#else
    nc_catch(nc_put_vara_double(output->ncid, output->varid[id], &start, &count, data))
#endif
    (void)time;
    (void)column;
    return;
}
