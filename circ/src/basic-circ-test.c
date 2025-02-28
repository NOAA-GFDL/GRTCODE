#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "argparse.h"
#include "circ.h"
#include "driver.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"

#include "circ1.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define NUM_MOLECULES 7
#define NUM_CFCS 3
#define NUM_CIAS 3
#define NUM_CIA_SPECIES 2


typedef struct Output
{
    char * variables[36];
    SpectralGrid_t const * lw_grid;
    SpectralGrid_t const * sw_grid;
} Output_t;


/*Copy data from one array to another.*/
void copy_data(double const * source, fp_t * destination, size_t size, fp_t const * constant)
{
    size_t i;
    for (i=0; i<size; ++i)
    {
        destination[i] = (fp_t)(source[i]);
    }
    if (constant == NULL)
    {
        return;
    }
    for (i=0; i<size; ++i)
    {
        destination[i] *= *constant;
    }
    return;
}


/*Linearly interpolate ppmv values in pressure space.*/
void pressure_interpolate(fp_t * ppmv, double const * abundance, int num_layers,
                          fp_t * layer_pressure, fp_t * level_pressure)
{
    fp_t const to_ppmv = 1.e6;
    ppmv[0] = (fp_t)abundance[0]*to_ppmv;
    ppmv[num_layers] = (fp_t)abundance[num_layers - 1]*to_ppmv;
    int i;
    for (i=1; i<num_layers; ++i)
    {
        ppmv[i] = (abundance[i - 1] + (abundance[i] - abundance[i - 1])*
                  (level_pressure[i] - layer_pressure[i - 1])/
                  (layer_pressure[i] - layer_pressure[i - 1]));
        ppmv[i] *= to_ppmv;
    }
    return;
}


/*Reserve memory and read in atmospheric data.*/
Atmosphere_t create_atmosphere(Parser_t * parser)
{
    /*Add/parse command line arguments.*/
    snprintf(parser->description, desclen,
             "Calculates the radiative fluxes for the NASA CIRC test cases.");
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
    parse_args(*parser);

    /*Determine the number of levels.*/
    Atmosphere_t atm;
    atm.num_levels = num_levels;
    atm.num_layers = atm.num_levels - 1;
    atm.num_columns = 1;
    atm.num_times = 1;

    /*Pressure.*/
    atm.level_pressure = (fp_t *)malloc(sizeof(*(atm.level_pressure))*atm.num_levels);
    copy_data(level_pressure, atm.level_pressure, atm.num_levels, NULL);
    atm.layer_pressure = (fp_t *)malloc(sizeof(*(atm.layer_pressure))*atm.num_layers);
    copy_data(layer_pressure, atm.layer_pressure, atm.num_layers, NULL);

    /*Temperature.*/
    atm.level_temperature = (fp_t *)malloc(sizeof(*(atm.level_temperature))*atm.num_levels);
    copy_data(level_temperature, atm.level_temperature, atm.num_levels, NULL);
    atm.layer_temperature = (fp_t *)malloc(sizeof(*(atm.layer_temperature))*atm.num_layers);
    copy_data(layer_temperature, atm.layer_temperature, atm.num_layers, NULL);
    atm.surface_temperature = (fp_t *)malloc(sizeof(*(atm.surface_temperature)));
    atm.surface_temperature[0] = (fp_t)surface_temperature;

    /*Solar zenith angle.*/
    atm.solar_zenith_angle = (fp_t *)malloc(sizeof(*(atm.solar_zenith_angle)));
    atm.solar_zenith_angle[0] = (fp_t)cos(2.*M_PI*solar_zenith_angle/360.);

    /*Solar irradiance.*/
    atm.total_solar_irradiance = (fp_t *)malloc(sizeof(*(atm.total_solar_irradiance)));
    atm.total_solar_irradiance[0] = toa_solar_irradiance/atm.solar_zenith_angle[0];

    /*Surface albedo.*/
    char buffer[valuelen];
    if (get_argument(*parser, "-a", buffer))
    {
        atm.albedo_grid_size = 2;
        atm.albedo_grid = (fp_t *)malloc(sizeof(*(atm.albedo_grid))*atm.albedo_grid_size);
        atm.albedo_grid[0] = -1.;
        atm.albedo_grid[1] = 0.;
        atm.surface_albedo =(fp_t *)malloc(sizeof(*(atm.surface_albedo))*atm.albedo_grid_size);
        atm.surface_albedo[0] = (fp_t)(atof(buffer));
        atm.surface_albedo[1] = atm.surface_albedo[0];
    }
    else
    {
        atm.albedo_grid_size = num_wavenumbers;
        atm.albedo_grid = (fp_t *)malloc(sizeof(*(atm.albedo_grid))*atm.albedo_grid_size);
        copy_data(wavenumber, atm.albedo_grid, atm.albedo_grid_size, NULL);
        atm.surface_albedo = (fp_t *)malloc(sizeof(*(atm.surface_albedo))*atm.albedo_grid_size);
        copy_data(surface_albedo, atm.surface_albedo, atm.albedo_grid_size, NULL);
    }

    /*Surface emissivity. CIRC cases assume this is 1.*/
    atm.emissivity_grid_size = 2;
    atm.emissivity_grid = (fp_t *)malloc(sizeof(*(atm.emissivity_grid))*atm.emissivity_grid_size);
    atm.emissivity_grid[0] = -1.;
    atm.emissivity_grid[1] = 0.;
    atm.surface_emissivity = (fp_t *)malloc(sizeof(*(atm.surface_emissivity))*atm.emissivity_grid_size);
    atm.surface_emissivity[0] = 1.;
    atm.surface_emissivity[1] = atm.surface_emissivity[0];

    /*Layer thickness.*/
    atm.layer_thickness = (fp_t *)malloc(sizeof(*(atm.layer_thickness))*atm.num_layers);
    int i;
    for (i=0; i<atm.num_layers; ++i)
    {
        fp_t const gas_constant = 8.314462; /*[J mol-1 K-1].*/
        fp_t const gravity = 9.81; /*[m s-2].*/
        fp_t const kg_per_g = .001; /*[kg g-1].*/
        fp_t const molar_mass = 28.9647; /*[g mol-1].*/
        atm.layer_thickness[i] = (fabs(log(atm.level_pressure[i]) - log(atm.level_pressure[i + 1]))*
                                 atm.layer_temperature[i]*gas_constant)/(molar_mass*kg_per_g*
                                 gravity);
    }

    /*Molecular abundances.*/
    struct MoleculeMeta
    {
        int id;
        char * flag;
        double * data;
    };
    int const num_molecules = NUM_MOLECULES;
    struct MoleculeMeta molecules[NUM_MOLECULES] = {
        {CH4, "-CH4", CH4_abundance},
        {CO, "-CO", CO_abundance},
        {CO2, "-CO2", CO2_abundance},
        {H2O, "-H2O", H2O_abundance},
        {N2O, "-N2O", N2O_abundance},
        {O2, "-O2", O2_abundance},
        {O3, "-O3", O3_abundance}
    };
    atm.molecules = (int *)malloc(sizeof(*(atm.molecules))*num_molecules);
    atm.num_molecules = 0;
    atm.ppmv = (fp_t **)malloc(sizeof(*(atm.ppmv))*num_molecules);
    for (i=0; i<num_molecules; ++i)
    {
        if (get_argument(*parser, molecules[i].flag, NULL))
        {
            atm.molecules[atm.num_molecules] = molecules[i].id;
            atm.ppmv[atm.num_molecules] = (fp_t *)malloc(sizeof(fp_t)*atm.num_levels);
            pressure_interpolate(atm.ppmv[atm.num_molecules], molecules[i].data,
                                 atm.num_layers, atm.layer_pressure, atm.level_pressure);
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
    int const num_cfcs = NUM_CFCS;
    struct MoleculeMeta cfcs[NUM_CFCS] = {
        {CFC11, "-CFC-11", CFC11_abundance},
        {CFC12, "-CFC-12", CFC12_abundance},
        {CCl4, "-CCl4", CCl4_abundance}
    };
    atm.cfc = (Cfc_t *)malloc(sizeof(*(atm.cfc))*num_cfcs);
    atm.num_cfcs = 0;
    atm.cfc_ppmv = (fp_t **)malloc(sizeof(*(atm.cfc_ppmv))*num_cfcs);
    for (i=0; i<num_cfcs; ++i)
    {
        if (get_argument(*parser, cfcs[i].flag, atm.cfc[atm.num_cfcs].path))
        {
            atm.cfc[atm.num_cfcs].id = cfcs[i].id;
            atm.cfc_ppmv[atm.num_cfcs] = (fp_t *)malloc(sizeof(fp_t)*atm.num_levels);
            pressure_interpolate(atm.cfc_ppmv[atm.num_cfcs], cfcs[i].data,
                                 atm.num_layers, atm.layer_pressure, atm.level_pressure);
            atm.num_cfcs++;
        }
    }

    /*Collision-induced absorption (CIA) abundances.*/
    int const num_cias = NUM_CIAS;
    int const num_cia_species = NUM_CIA_SPECIES;
    struct CiaMeta
    {
        int species1;
        int species2;
        char *flag;
    };
    struct CiaMeta cias[NUM_CIAS] = {
        {CIA_N2, CIA_N2, "-N2-N2"},
        {CIA_O2, CIA_N2, "-O2-N2"},
        {CIA_O2, CIA_O2, "-O2-O2"}
    };
    atm.cia = (Cia_t *)malloc(sizeof(*(atm.cia))*num_cias);
    atm.num_cias = 0;
    atm.cia_species = (int *)malloc(sizeof(*(atm.cia_species))*num_cia_species);
    atm.num_cia_species = 0;
    atm.cia_ppmv = (fp_t **)malloc(sizeof(*(atm.cia_ppmv))*num_cia_species);
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
                    atm.cia_ppmv[atm.num_cia_species] = (fp_t *)malloc(sizeof(fp_t)*atm.num_levels);
                    fp_t * ppmv = atm.cia_ppmv[atm.num_cia_species];
                    if (atm.cia_species[atm.num_cia_species] == CIA_N2)
                    {
                        for (k=0; k<atm.num_levels; ++k)
                        {
                            ppmv[k] = 0.781*1.e6;
                        }
                    }
                    else if (atm.cia_species[atm.num_cia_species] == CIA_O2)
                    {
                        pressure_interpolate(ppmv, O2_abundance, atm.num_layers,
                                             atm.layer_pressure, atm.level_pressure);
                    }
                    atm.num_cia_species++;
                }
            }
            atm.num_cias++;
        }
    }

    /*Aerosols.*/
    atm.clean = get_argument(*parser, "-clean", NULL);
    if (!atm.clean)
    {
        fp_t alpha = angstrom_exponent_value;
        atm.aerosol_grid_size = 50000;
        atm.aerosol_grid = (fp_t *)malloc(sizeof(*(atm.aerosol_grid))*atm.aerosol_grid_size);
        for (i=0; i<atm.aerosol_grid_size; ++i)
        {
            atm.aerosol_grid[i] = (fp_t)(i + 1);
        }
        double * optics = aerosol_optical_depth_at_1_micron;
        atm.aerosol_optical_depth = (fp_t *)malloc(sizeof(*(atm.aerosol_optical_depth))*
                                                   atm.num_layers*atm.aerosol_grid_size);
        fp_t const cmtomicron = 10000.;
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                fp_t lambda = cmtomicron/(atm.aerosol_grid[j]);
                atm.aerosol_optical_depth[i*atm.aerosol_grid_size + j] = optics[i]*pow(lambda, -1.*alpha);
            }
        }
        optics = aerosol_single_scatter_albedo;
        atm.aerosol_single_scatter_albedo = (fp_t *)malloc(sizeof(*(atm.aerosol_single_scatter_albedo))*
                                                           atm.num_layers*atm.aerosol_grid_size);
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                atm.aerosol_single_scatter_albedo[i*atm.aerosol_grid_size + j] = optics[i];
            }
        }
        optics = aerosol_asymmetry_factor;
        atm.aerosol_asymmetry_factor = (fp_t *)malloc(sizeof(*(atm.aerosol_asymmetry_factor))*
                                                      atm.num_layers*atm.aerosol_grid_size);
        for (i=0; i<atm.num_layers; ++i)
        {
            int j;
            for (j=0; j<atm.aerosol_grid_size; ++j)
            {
                atm.aerosol_asymmetry_factor[i*atm.aerosol_grid_size + j] = optics[i];
            }
        }
    }

    /*Clouds.*/
    atm.clear = get_argument(*parser, "-clear", NULL);
    if (!atm.clear)
    {
        atm.cloud_fraction = (fp_t *)malloc(sizeof(*(atm.cloud_fraction))*atm.num_layers);
        copy_data(cloud_fraction, atm.cloud_fraction, atm.num_layers, NULL);
        atm.ice_water_content = (fp_t *)malloc(sizeof(*(atm.ice_water_content))*atm.num_layers);
        atm.liquid_water_content = (fp_t *)malloc(sizeof(*(atm.liquid_water_content))*atm.num_layers);
        for (i=0; i<atm.num_layers; ++i)
        {
            atm.ice_water_content[i] = ice_water_path[i]/atm.layer_thickness[i];
            atm.liquid_water_content[i] = liquid_water_path[i]/atm.layer_thickness[i];
        }
        atm.liquid_water_droplet_radius = (fp_t *)malloc(sizeof(*(atm.liquid_water_droplet_radius))*atm.num_layers);
        copy_data(liquid_water_effective_particle_size, atm.liquid_water_droplet_radius,
                  atm.num_layers, NULL);
    }
    return atm;
}


/*Free memory for atmosphere.*/
void destroy_atmosphere(Atmosphere_t * const atm)
{
    if (!atm->clean)
    {
        free(atm->aerosol_grid);
        free(atm->aerosol_optical_depth);
        free(atm->aerosol_single_scatter_albedo);
        free(atm->aerosol_asymmetry_factor);
    }
    free(atm->albedo_grid);
    free(atm->cfc);
    int i;
    for (i=0; i<atm->num_cfcs; ++i)
    {
        free(atm->cfc_ppmv[i]);
    }
    free(atm->cfc_ppmv);
    free(atm->cia);
    for (i=0; i<atm->num_cia_species; ++i)
    {
        free(atm->cia_ppmv[i]);
    }
    free(atm->cia_ppmv);
    free(atm->cia_species);
    free(atm->emissivity_grid);
    free(atm->layer_pressure);
    free(atm->layer_temperature);
    free(atm->layer_thickness);
    free(atm->level_pressure);
    free(atm->level_temperature);
    free(atm->molecules);
    for (i=0; i<atm->num_molecules; ++i)
    {
        free(atm->ppmv[i]);
    }
    free(atm->ppmv);
    free(atm->solar_zenith_angle);
    free(atm->surface_albedo);
    free(atm->surface_emissivity);
    free(atm->surface_temperature);
    free(atm->total_solar_irradiance);
    if (!atm->clear)
    {
        free(atm->cloud_fraction);
        free(atm->ice_water_content);
        free(atm->liquid_water_content);
        free(atm->liquid_water_droplet_radius);
    }
    return;
}


/*Create an output file and write metadata.*/
void create_flux_file(Output_t ** output,
                      char const * const filepath,
                      Atmosphere_t const * const atm,
                      SpectralGrid_t const * const lw_grid,
                      SpectralGrid_t const * const sw_grid
)
{
    Output_t * file = (Output_t *)malloc(sizeof(*file));
    file->lw_grid = lw_grid;
    file->sw_grid = sw_grid;
    *output = file;
    return;
}


/*Close output file.*/
void close_flux_file(Output_t * const output)
{
    free(output);
    return;
}


/*Write fluxes to the output file.*/
void write_output(Output_t * output, unsigned int id, fp_t const * data, int time, int column)
{
    char * buffer;
    SpectralGrid_t const * grid;
    fp_t lblrtm_reference;
    fp_t mean_reference;
    if (id == (LONGWAVE ^ UPWARD ^ TOP ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RLUTCSAF";
        grid = output->lw_grid;
        lblrtm_reference = 304.27;
        mean_reference = 301.7;
    }
    else if (id == (LONGWAVE ^ UPWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RLUSCSAF";
        grid = output->lw_grid;
        lblrtm_reference = 445.12;
        mean_reference = 0.;
    }
    else if (id == (LONGWAVE ^ DOWNWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RLDSCSAF";
        grid = output->lw_grid;
        lblrtm_reference = 288.2;
        mean_reference = 289.7;
    }
    else if (id == (SHORTWAVE ^ UPWARD ^ TOP ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RSUTCSAF";
        grid = output->sw_grid;
        lblrtm_reference = 175.0;
        mean_reference = 169.8;
    }
    else if (id == (SHORTWAVE ^ UPWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RSUSCSAF";
        grid = output->sw_grid;
        lblrtm_reference = 137.40;
        mean_reference = 0;
    }
    else if (id == (SHORTWAVE ^ DOWNWARD ^ TOP ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RSDTCSAF";
        grid = output->sw_grid;
        lblrtm_reference = 912.79;
        mean_reference = 0;
    }
    else if (id == (SHORTWAVE ^ DOWNWARD ^ SURFACE ^ CLEARSKY ^ AEROSOLFREE))
    {
        buffer = "RSDSCSAF";
        grid = output->sw_grid;
        lblrtm_reference = 701.2;
        mean_reference = 705.9;
    }
    else
    {
        return;
    }
    int i;
    fp_t integrated_flux = 0.;
    for (i=0; i<grid->n - 1; ++i)
    {
        integrated_flux += 0.5*(data[i] + data[i + 1])*grid->dw;
    }
    fprintf(stdout, "%s: %e\t\t%e\t\t%e\n", buffer, lblrtm_reference, mean_reference,
            integrated_flux);
    (void)time;
    return;
}
