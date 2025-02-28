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
#include "rfmip-irf.h"
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
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


struct Output
{
    int dimid[256];
    int ncid;
    int varid[256];
    SpectralGrid_t const * lw_grid;
    SpectralGrid_t const * sw_grid;
    int integrated; /*Are the output values integrated (1) or spectrally resolved (0).*/
};


/*Reserve memory and read in atmospheric data.*/
Atmosphere_t create_atmosphere(Parser_t *parser)
{
    /*Add/parse command line arguments.*/
    snprintf(parser->description, desclen,
             "Calculates the radiative fluxes for the CMIP6 RFMIP-IRF scenarios.");
    add_argument(parser, "input_file", NULL, "Input data file.", NULL);
    add_argument(parser, "experiment", NULL, "Experiment number.", NULL);
    int one = 1;
    add_argument(parser, "-CCl4", NULL, "CSV file with CCl4 cross sections.", &one);
    add_argument(parser, "-C2F6", NULL, "CSV file with C2F6 cross sections.", &one);
    add_argument(parser, "-CF4", NULL, "CSV file with CF4 cross sections.", &one);
    add_argument(parser, "-CFC-11", NULL, "CSV file with CFC-11 cross sections.", &one);
    add_argument(parser, "-CFC-11-eq", NULL, "CSV file with CFC-11 cross sections.", &one);
    add_argument(parser, "-CFC-12", NULL, "CSV file with CFC-12 cross sections.", &one);
    add_argument(parser, "-CFC-12-eq", NULL, "CSV file with CFC-12 cross sections.", &one);
    add_argument(parser, "-CFC-113", NULL, "CSV file with CFC-113 cross sections.", &one);
    add_argument(parser, "-CFC-114", NULL, "CSV file with CFC-114 cross sections.", &one);
    add_argument(parser, "-CFC-115", NULL, "CSV file with CFC-115 cross sections.", &one);
    add_argument(parser, "-CH2Cl2", NULL, "CSV file with CH2Cl2 cross sections.", &one);
    add_argument(parser, "-CH4", NULL, "Include CH4.", NULL);
    add_argument(parser, "-CO", NULL, "Include CO.", NULL);
    add_argument(parser, "-CO2", NULL, "Include CO2.", NULL);
    add_argument(parser, "-H2O", NULL, "Include H2O.", NULL);
    add_argument(parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(parser, "-HCFC-22", NULL, "CSV file with HCFC-22 cross sections.", &one);
    add_argument(parser, "-HCFC-141b", NULL, "CSV file with HCFC-141b cross sections.", &one);
    add_argument(parser, "-HCFC-142b", NULL, "CSV file with HCFC-142b cross sections.", &one);
    add_argument(parser, "-HFC-23", NULL, "CSV file with HFC-23 cross sections.", &one);
    add_argument(parser, "-HFC-125", NULL, "CSV file with HFC-125 cross sections.", &one);
    add_argument(parser, "-HFC-134a", NULL, "CSV file with HFC-134a cross sections.", &one);
    add_argument(parser, "-HFC-134a-eq", NULL, "CSV file with HFC-134a cross sections.", &one);
    add_argument(parser, "-HFC-143a", NULL, "CSV file with HFC-143a cross sections.", &one);
    add_argument(parser, "-HFC-152a", NULL, "CSV file with HFC-152a cross sections.", &one);
    add_argument(parser, "-HFC-227ea", NULL, "CSV file with HFC-227ea cross sections.", &one);
    add_argument(parser, "-HFC-245fa", NULL, "CSV file with HFC-245fa cross sections.", &one);
    add_argument(parser, "-N2-N2", NULL, "CSV file with N2-N2 collison cross sections", &one);
    add_argument(parser, "-N2O", NULL, "Include N2O.", NULL);
    add_argument(parser, "-NF3", NULL, "CSV file with NF3 cross sections.", &one);
    add_argument(parser, "-O2", NULL, "Include O2.", NULL);
    add_argument(parser, "-O2-N2", NULL, "CSV file with O2-N2 collison cross sections", &one);
    add_argument(parser, "-O2-O2", NULL, "CSV file with O2-O2 collison cross sections", &one);
    add_argument(parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(parser, "-SF6", NULL, "CSV file with SF6 cross sections.", &one);
    add_argument(parser, "-x", "--column-lower-bound", "Starting column index.", &one);
    add_argument(parser, "-X", "--column-upper-bound", "Ending column index.", &one);
    add_argument(parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(*parser);

    /*Open the input file.*/
    char buffer[valuelen];
    get_argument(*parser, "input_file", buffer);
    int ncid;
    nc_catch(nc_open(buffer, NC_NOWRITE, &ncid));

    /*Determine the experiment to run.*/
    get_argument(*parser, "experiment", buffer);
    int experiment = atoi(buffer);

    /*Determine the number of columns.*/
    Atmosphere_t atm;
    int x = get_argument(*parser, "-x", buffer) ? atoi(buffer) : 0;
    int X;
    if (get_argument(*parser, "-X", buffer))
    {
        X = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "site", &dimid));
        char name[128];
        memset(name, 0, 128*sizeof(char));
        size_t num_columns;
        nc_catch(nc_inq_dim(ncid, dimid, name, &num_columns));
        X = (int)num_columns - 1;
    }
    atm.x = x;
    atm.X = X;
    atm.num_columns = X - x + 1;

    /*Determine the number of levels.*/
    int z = get_argument(*parser, "-z", buffer) ? atoi(buffer) : 0;
    int Z;
    if (get_argument(*parser, "-Z", buffer))
    {
        Z = atoi(buffer);
    }
    else
    {
        int dimid;
        nc_catch(nc_inq_dimid(ncid, "level", &dimid));
        char name[128];
        memset(name, 0, 128*sizeof(char));
        size_t num_levels;
        nc_catch(nc_inq_dim(ncid, dimid, name, &num_levels));
        Z = (int)num_levels - 1;
    }
    atm.num_levels = Z - z + 1;
    atm.num_layers = atm.num_levels - 1;
    atm.num_times = 1;

    /*Pressure.*/
    alloc(atm.level_pressure, atm.num_levels*atm.num_columns, fp_t *);
    int varid;
    nc_catch(nc_inq_varid(ncid, "pres_level", &varid));
    size_t start[3] = {x, z, 0};
    size_t count[3] = {atm.num_columns, atm.num_levels, 1};
    get_var(ncid, varid, start, count, atm.level_pressure);
    fp_t const Patomb = 0.01;
    int i;
    for (i=0; i<atm.num_columns*atm.num_levels; ++i)
    {
        atm.level_pressure[i] *= Patomb;
    }
    alloc(atm.layer_pressure, atm.num_layers*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "pres_layer", &varid));
    start[0] = x; start[1] = z; start[2] = 0;
    count[0] = atm.num_columns; count[1] = atm.num_layers; count[2] = 1;
    get_var(ncid, varid, start, count, atm.layer_pressure);
    for (i=0; i<atm.num_columns*atm.num_layers; ++i)
    {
        atm.layer_pressure[i] *= Patomb;
    }

    /*Temperature.*/
    alloc(atm.level_temperature, atm.num_levels*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "temp_level", &varid));
    start[0] = experiment; start[1] = x; start[2] = z;
    count[0] = 1; count[1] = atm.num_columns; count[2] = atm.num_levels;
    get_var(ncid, varid, start, count, atm.level_temperature);
    alloc(atm.layer_temperature, atm.num_layers*atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "temp_layer", &varid));
    start[0] = experiment; start[1] = x; start[2] = z;
    count[0] = 1; count[1] = atm.num_columns; count[2] = atm.num_layers;
    get_var(ncid, varid, start, count, atm.layer_temperature);
    alloc(atm.surface_temperature, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_temperature", &varid));
    start[0] = experiment; start[1] = x; start[2] = 0;
    count[0] = 1; count[1] = atm.num_columns; count[2] = 1;
    get_var(ncid, varid, start, count, atm.surface_temperature);

    /*Solar zenith angle.*/
    alloc(atm.solar_zenith_angle, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "solar_zenith_angle", &varid));
    start[0] = x; start[1] = 0; start[2] = 0;
    count[0] = atm.num_columns; count[1] = 1; count[2] = 1;
    get_var(ncid, varid, start, count, atm.solar_zenith_angle);
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.solar_zenith_angle[i] = cos(2.*M_PI*atm.solar_zenith_angle[i]/360.);
    }

    /*Solar irradiance.*/
    alloc(atm.total_solar_irradiance, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "total_solar_irradiance", &varid));
    start[0] = x; start[1] = 0; start[2] = 0;
    count[0] = atm.num_columns; count[1] = 1; count[2] = 1;
    get_var(ncid, varid, start, count, atm.total_solar_irradiance);

    /*Surface albedo.*/
    atm.albedo_grid_size = 2;
    alloc(atm.albedo_grid, atm.albedo_grid_size, fp_t *);
    atm.albedo_grid[0] = -1.;
    atm.albedo_grid[1] = 0.;
    alloc(atm.surface_albedo, atm.num_columns*atm.albedo_grid_size, fp_t *);
    fp_t *albedo;
    alloc(albedo, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_albedo", &varid));
    start[0] = x; start[1] = 0; start[2] = 0;
    count[0] = atm.num_columns; count[1] = 1; count[2] = 1;
    get_var(ncid, varid, start, count, albedo);
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.surface_albedo[i*atm.albedo_grid_size] = albedo[i];
        atm.surface_albedo[i*atm.albedo_grid_size+1] = albedo[i];
    }
    free(albedo);

    /*Surface emissivity.*/
    atm.emissivity_grid_size = 2;
    alloc(atm.emissivity_grid, atm.emissivity_grid_size, fp_t *);
    atm.emissivity_grid[0] = -1.;
    atm.emissivity_grid[1] = 0.;
    alloc(atm.surface_emissivity, atm.num_columns*atm.emissivity_grid_size, fp_t *);
    fp_t *emissivity;
    alloc(emissivity, atm.num_columns, fp_t *);
    nc_catch(nc_inq_varid(ncid, "surface_emissivity", &varid));
    start[0] = x; start[1] = 0; start[2] = 0;
    count[0] = atm.num_columns; count[1] = 1; count[2] = 1;
    get_var(ncid, varid, start, count, emissivity);
    for (i=0; i<atm.num_columns; ++i)
    {
        atm.surface_emissivity[i*atm.emissivity_grid_size] = emissivity[i];
        atm.surface_emissivity[i*atm.emissivity_grid_size+1] = emissivity[i];
    }
    free(emissivity);

    /*Molecular abundances.*/
    fp_t const to_ppmv = 1.e6;
    struct MoleculeMeta
    {
        int id;
        char *flag;
        char *name;
        int global_mean;
    };
    int const num_molecules = 7;
    struct MoleculeMeta molecules[32] = {
        {CH4, "-CH4", "methane_GM", 1},
        {CO, "-CO", "carbon_monoxide_GM", 1},
        {CO2, "-CO2", "carbon_dioxide_GM", 1},
        {H2O, "-H2O", "water_vapor", 0},
        {N2O, "-N2O", "nitrous_oxide_GM", 1},
        {O2, "-O2", "oxygen_GM", 1},
        {O3, "-O3", "ozone", 0}
    };
    alloc(atm.molecules, num_molecules, int *);
    atm.num_molecules = 0;
    alloc(atm.ppmv, num_molecules, fp_t **);
    fp_t *abundance;
    alloc(abundance, atm.num_columns*atm.num_layers, fp_t *);
    for (i=0; i<num_molecules; ++i)
    {
        if (get_argument(*parser, molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = molecules[i].id;
            alloc(atm.ppmv[atm.num_molecules], atm.num_levels*atm.num_columns, fp_t *);
            fp_t *ppmv = atm.ppmv[atm.num_molecules];
            if (!molecules[i].global_mean)
            {
                nc_catch(nc_inq_varid(ncid, molecules[i].name, &varid));
                start[0] = experiment; start[1] = x; start[2] = z;
                count[0] = 1; count[1] = atm.num_columns; count[2] = atm.num_layers;
                get_var(ncid, varid, start, count, abundance);
                int j;
                for (j=0; j<atm.num_columns; ++j)
                {
                    int o = j*atm.num_layers;
                    ppmv[j*atm.num_levels] = abundance[o]*to_ppmv;
                    ppmv[(j+1)*atm.num_levels-1] = abundance[o+atm.num_layers-1]*to_ppmv;
                    int k;
                    for (k=1; k<atm.num_layers; ++k)
                    {
                        ppmv[j*atm.num_levels+k] = to_ppmv*(abundance[o+k-1] + (abundance[o+k] - abundance[o+k-1])*
                                                   (atm.level_pressure[j*atm.num_levels+k] - atm.layer_pressure[o+k-1])/
                                                   (atm.layer_pressure[o+k] - atm.layer_pressure[o+k-1]));
                    }
                }
            }
            else
            {
                fp_t gm = 0.;
                nc_catch(nc_inq_varid(ncid, molecules[i].name, &varid));
                start[0] = experiment; start[1] = 0; start[2] = 0;
                count[0] = 1; count[1] = 1; count[2] = 1;
                get_var(ncid, varid, start, count, &gm);
                char units[32];
                memset(units, 0, 32*sizeof(char));
                nc_catch(nc_get_att_text(ncid, varid, "units", units));
                gm *= atof(units)*to_ppmv;
                int j;
                for (j=0; j<atm.num_columns*atm.num_levels; ++j)
                {
                    ppmv[j] = gm;
                }
            }
            atm.num_molecules++;
        }
    }
    free(abundance);

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
    int const num_cfcs = 24;
    struct MoleculeMeta cfcs[32] = {
        {CCl4, "-CCl4", "carbon_tetrachloride_GM", 1},
        {C2F6, "-C2F6", "c2f6_GM", 1},
        {CF4, "-CF4", "cf4_GM", 1},
        {CFC11, "-CFC-11", "cfc11_GM", 1},
        {CFC11, "-CFC-11-eq", "cfc11eq_GM", 1},
        {CFC12, "-CFC-12", "cfc12_GM", 1},
        {CFC12, "-CFC-12-eq", "cfc12eq_GM", 1},
        {CFC113, "-CFC-113", "cfc113_GM", 1},
        {CFC114, "-CFC-114", "cfc114_GM", 1},
        {CFC115, "-CFC-115", "cfc115_GM", 1},
        {CH2Cl2, "-CH2Cl2", "ch2cl2_GM", 1},
        {HCFC22, "-HCFC-22", "hcfc22_GM", 1},
        {HCFC141b, "-HCFC-141b", "hcfc141b_GM", 1},
        {HCFC142b, "-HCFC-142b", "hcfc142b_GM", 1},
        {HFC23, "-HFC-23", "hfc23_GM", 1},
        {HFC125, "-HFC-125", "hfc125_GM", 1},
        {HFC134a, "-HFC-134a", "hfc134a_GM", 1},
        {HFC134a, "-HFC-134a-eq", "hfc134aeq_GM", 1},
        {HFC143a, "-HFC-143a", "hfc143a_GM", 1},
        {HFC152a, "-HFC-152a", "hfc152a_GM", 1},
        {HFC227ea, "-HFC-227ea", "hfc227ea_GM", 1},
        {HFC245fa, "-HFC-245fa", "hfc245fa_GM", 1},
        {NF3, "-NF3", "nf3_GM", 1},
        {SF6, "-SF6", "sf6_GM", 1}
    };
    alloc(atm.cfc, num_cfcs, Cfc_t *);
    atm.num_cfcs = 0;
    alloc(atm.cfc_ppmv, num_cfcs, fp_t **);
    for (i=0; i<num_cfcs; ++i)
    {
        if (get_argument(*parser, cfcs[i].flag, atm.cfc[atm.num_cfcs].path))
        {
            atm.cfc[atm.num_cfcs].id = cfcs[i].id;
            alloc(atm.cfc_ppmv[atm.num_cfcs], atm.num_levels*atm.num_columns, fp_t *);
            fp_t *ppmv = atm.cfc_ppmv[atm.num_cfcs];
            fp_t gm;
            nc_catch(nc_inq_varid(ncid, cfcs[i].name, &varid));
            start[0] = experiment; start[1] = 0; start[2] = 0;
            count[0] = 1; count[1] = 1; count[2] = 1;
            get_var(ncid, varid, start, count, &gm);
            char units[32];
            memset(units, 0, 32*sizeof(char));
            nc_catch(nc_get_att_text(ncid, varid, "units", units));
            gm *= atof(units)*to_ppmv;
            int j;
            for (j=0; j<atm.num_columns*atm.num_levels; ++j)
            {
                ppmv[j] = gm;
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
    struct CiaMeta cias[8] = {
        {CIA_N2, CIA_N2, "-N2-N2"},
        {CIA_O2, CIA_N2, "-O2-N2"},
        {CIA_O2, CIA_O2, "-O2-O2"}
    };
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
                    char *name;
                    if (atm.cia_species[atm.num_cia_species] == CIA_N2)
                    {
                        name = "nitrogen_GM";
                    }
                    else if (atm.cia_species[atm.num_cia_species] == CIA_O2)
                    {
                        name = "oxygen_GM";
                    }
                    fp_t gm;
                    nc_catch(nc_inq_varid(ncid, name, &varid));
                    start[0] = experiment; start[1] = 0; start[2] = 0;
                    count[0] = 1; count[1] = 1; count[2] = 1;
                    get_var(ncid, varid, start, count, &gm);
                    char units[32];
                    memset(units, 0, 32*sizeof(char));
                    nc_catch(nc_get_att_text(ncid, varid, "units", units));
                    gm *= atof(units)*to_ppmv;
                    int j;
                    for (j=0; j<atm.num_columns*atm.num_levels; ++j)
                    {
                        ppmv[j] = gm;
                    }
                    atm.num_cia_species++;
                }
            }
            atm.num_cias++;
        }
    }

    /*Aerosols.*/
    atm.clean = 1;

    /*Clouds.*/
    atm.clear = 1;

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
static void add_variable(Output_t * const o, /*Output object.*/
                         Variables_t const index, /*Variable index.*/
                         char const * const name, /*Variable name.*/
                         char const * const standard_name, /*Variable standard name.*/
                         char const * const units, /*Units.*/
                         fp_t const * const fill_value, /*Fill value.*/
                         int integrated,
                         Variables_t wavenumber
)
{
#ifdef SINGLE_PRESCISION
    nc_type type = NC_FLOAT;
#else
    nc_type type = NC_DOUBLE;
#endif
    int varid;
    if (integrated)
    {
        nc_catch(nc_def_var(o->ncid, name, type, 1, &(o->dimid[COLUMN]), &varid));
    }
    else
    {
        int dimids[2] = {o->dimid[COLUMN], o->dimid[wavenumber]};
        nc_catch(nc_def_var(o->ncid, name, type, 2, dimids, &varid));
    }
    nc_catch(nc_put_att_text(o->ncid, varid, "units", strlen(units), units));
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
void create_flux_file(Output_t ** output, char const * const filepath,
                      Atmosphere_t const * const atm, SpectralGrid_t const * const lw_grid,
                      SpectralGrid_t const * const sw_grid, int user_level,
                      int integrated)
{
    Output_t * file = (Output_t *)malloc(sizeof(*file));
    int i;
    for (i=0; i<256; ++i)
    {
        file->dimid[i] = -1;
        file->varid[i] = -1;
    }
    nc_catch(nc_create(filepath, NC_NETCDF4, &(file->ncid)));
    nc_catch(nc_put_att_int(file->ncid, NC_GLOBAL, "x_start", NC_INT, 1, &(atm->x)));
    nc_catch(nc_put_att_int(file->ncid, NC_GLOBAL, "x_stop", NC_INT, 1, &(atm->X)));
    nc_catch(nc_def_dim(file->ncid, "level", atm->num_levels, &(file->dimid[LEVEL])));
    nc_catch(nc_def_var(file->ncid, "level", NC_DOUBLE, 1, &(file->dimid[LEVEL]),
                        &(file->varid[LEVEL])));
    nc_catch(nc_put_att_text(file->ncid, file->varid[LEVEL], "axis", 2, "z"));
    nc_catch(nc_def_dim(file->ncid, "column", atm->num_columns, &(file->dimid[COLUMN])));
    nc_catch(nc_def_var(file->ncid, "column", NC_DOUBLE, 1, &(file->dimid[COLUMN]),
                        &(file->varid[COLUMN])));
    nc_catch(nc_put_att_text(file->ncid, file->varid[COLUMN], "axis", 2, "x"));
    if (!integrated)
    {
        nc_catch(nc_def_dim(file->ncid, "lw_wavenumber", lw_grid->n, &(file->dimid[LW_WAVENUMBER])));
        nc_catch(nc_def_var(file->ncid, "lw_wavenumber", NC_DOUBLE, 1, &(file->dimid[LW_WAVENUMBER]),
                            &(file->varid[LW_WAVENUMBER])));
        nc_catch(nc_put_att_text(file->ncid, file->varid[LW_WAVENUMBER], "units", 5, "cm-1"));
        nc_catch(nc_def_dim(file->ncid, "sw_wavenumber", sw_grid->n, &(file->dimid[SW_WAVENUMBER])));
        nc_catch(nc_def_var(file->ncid, "sw_wavenumber", NC_DOUBLE, 1, &(file->dimid[SW_WAVENUMBER]),
                            &(file->varid[SW_WAVENUMBER])));
        nc_catch(nc_put_att_text(file->ncid, file->varid[LW_WAVENUMBER], "units", 5, "cm-1"));
    }

    add_variable(file, RLUTCSAF, "rlutcsaf", "upwelling_toa_longwave_flux_in_air", "W m-2",
                 NULL, integrated, LW_WAVENUMBER);
    add_variable(file, RLUSCSAF, "rluscsaf", "upwelling_surface_longwave_flux_in_air", "W m-2",
                 NULL, integrated, LW_WAVENUMBER);
    add_variable(file, RLDSCSAF, "rldscsaf", "downwelling_surface_longwave_flux_in_air", "W m-2",
                 NULL, integrated, LW_WAVENUMBER);

    add_variable(file, RLDCSAF_USER_LEVEL, "rldcsaf_level", "downwelling_longwave_flux_in_air",
                 "W m-2", NULL, integrated, LW_WAVENUMBER);
    nc_catch(nc_put_att_int(file->ncid, file->varid[RLDCSAF_USER_LEVEL], "level",
                            NC_INT, 1, &user_level));

    add_variable(file, RLUCSAF_USER_LEVEL, "rlucsaf_level", "upwelling_longwave_flux_in_air",
                 "W m-2", NULL, integrated, LW_WAVENUMBER);
    nc_catch(nc_put_att_int(file->ncid, file->varid[RLUCSAF_USER_LEVEL], "level",
                            NC_INT, 1, &user_level));

    fp_t const zero = 0;
    add_variable(file, RSUTCSAF, "rsutcsaf", "upwelling_toa_shortwave_flux_in_air", "W m-2",
                 &zero, integrated, SW_WAVENUMBER);
    add_variable(file, RSUSCSAF, "rsuscsaf", "upwelling_surface_shortwave_flux_in_air", "W m-2",
                 &zero, integrated, SW_WAVENUMBER);
    add_variable(file, RSDTCSAF, "rsdtcsaf", "downwelling_toa_shortwave_flux_in_air", "W m-2",
                 &zero, integrated, SW_WAVENUMBER);
    add_variable(file, RSDSCSAF, "rsdscsaf", "downwelling_surface_shortwave_flux_in_air", "W m-2",
                 &zero, integrated, SW_WAVENUMBER);

    add_variable(file, RSDCSAF_USER_LEVEL, "rsdcsaf_level", "downwelling_shortwave_flux_in_air",
                 "W m-2", &zero, integrated, SW_WAVENUMBER);
    nc_catch(nc_put_att_int(file->ncid, file->varid[RSDCSAF_USER_LEVEL], "level",
                            NC_INT, 1, &user_level));

    add_variable(file, RSUCSAF_USER_LEVEL, "rsucsaf_level", "upwelling_shortwave_flux_in_air",
                 "W m-2", NULL, integrated, SW_WAVENUMBER);
    nc_catch(nc_put_att_int(file->ncid, file->varid[RSUCSAF_USER_LEVEL], "level",
                            NC_INT, 1, &user_level));

    file->lw_grid = lw_grid;
    file->sw_grid = sw_grid;
    file->integrated = integrated;
    *output = file;
    return;
}


/*Close output file.*/
void close_flux_file(Output_t * const o)
{
    nc_catch(nc_close(o->ncid));
    return;
}


/*Write fluxes to the output file.*/
void write_output(Output_t * output, Variables_t id, fp_t const * data, int time, int column)
{
    (void)time;
    if (id < 0)
    {
        /*This variable is not in the output file, so just silently return.*/
        return;
    }
    if (id > 256)
    {
        fprintf(stderr, "[%s: %d] Error: bad id = %d.\n", __FILE__, __LINE__, id); \
        exit(EXIT_FAILURE); \
    }
    if (output->varid[id] == -1)
    {
        return;
    }

    /*Get the data size.*/
    size_t start[2] = {column, 0};
    size_t count[2] = {1, 1};
    SpectralGrid_t const * grid;
    if (!output->integrated)
    {
        if (is_longwave_flux(id))
        {
            count[1] = output->lw_grid->n;
        }
        else if (is_shortwave_flux(id))
        {
            count[1] = output->sw_grid->n;
        }
    }

    /*Write out the data.*/
#ifdef SINGLE_PRECISION
    nc_catch(nc_put_vara_float(output->ncid, output->varid[id], start, count, data));
#else
    nc_catch(nc_put_vara_double(output->ncid, output->varid[id], start, count, data));
#endif
    return;
}
