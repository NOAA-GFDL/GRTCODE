#ifndef DRIVER_H_
#define DRIVER_H_

#include "argparse.h"
#include "grtcode_utilities.h"


typedef enum Variables
{
    MISSING = -1,
    TIME = 0,
    LATITUDE,
    LONGITUDE,
    COLUMN,
    LAYER,
    LEVEL,
    LW_WAVENUMBER,
    SW_WAVENUMBER,
    OPTICAL_DEPTH,
    VERTICALLY_INTEGRATED_OPTICAL_DEPTH,
    H2O_VMR,
    O3_VMR,
    CH4_VMR,
    CO2_VMR,
    N2O_VMR,
    LEVEL_TEMPERATURE,
    LAYER_TEMPERATURE,
    SURFACE_TEMPERATURE,
    LEVEL_PRESSURE,
    LAYER_PRESSURE,
    RLD,
    RLDAF,
    RLDCS,
    RLDCSAF,
    RLDS,
    RLDSAF,
    RLDSCS,
    RLDSCSAF,
    RLU,
    RLUAF,
    RLUCS,
    RLUCSAF,
    RLUS,
    RLUSAF,
    RLUSCS,
    RLUSCSAF,
    RLUT,
    RLUTAF,
    RLUTCS,
    RLUTCSAF,
    RLD_USER_LEVEL,
    RLDAF_USER_LEVEL,
    RLDCS_USER_LEVEL,
    RLDCSAF_USER_LEVEL,
    RLU_USER_LEVEL,
    RLUAF_USER_LEVEL,
    RLUCS_USER_LEVEL,
    RLUCSAF_USER_LEVEL,
    RSD,
    RSDAF,
    RSDCS,
    RSDCSAF,
    RSDS,
    RSDSAF,
    RSDSCS,
    RSDSCSAF,
    RSDT,
    RSDTAF,
    RSDTCS,
    RSDTCSAF,
    RSU,
    RSUAF,
    RSUCS,
    RSUCSAF,
    RSUS,
    RSUSAF,
    RSUSCS,
    RSUSCSAF,
    RSUT,
    RSUTAF,
    RSUTCS,
    RSUTCSAF,
    RSD_USER_LEVEL,
    RSDAF_USER_LEVEL,
    RSDCS_USER_LEVEL,
    RSDCSAF_USER_LEVEL,
    RSU_USER_LEVEL,
    RSUAF_USER_LEVEL,
    RSUCS_USER_LEVEL,
    RSUCSAF_USER_LEVEL,
    NUM_VARS
} Variables_t;


typedef struct Cfc
{
    int id; /*GRTCODE CFC id.*/
    char path[valuelen]; /*Path to CFC cross-section file.*/
} Cfc_t;


typedef struct Cia
{
    int id[2]; /*GRTCODE CIA id of species.*/
    char path[valuelen]; /*Path to CIA cross-section file.*/
} Cia_t;


int is_longwave_flux(Variables_t id);


int is_shortwave_flux(Variables_t id);


typedef struct Atmosphere
{
    fp_t *aerosol_grid; /*Wavenumber [cm-1] grid for the aerosol optical properties (wavenumber).*/
    int aerosol_grid_size; /*Length of aerosol_grid array.*/
    fp_t *aerosol_optical_depth; /*Aerosol optical depth (time, column, layer, wavenumber).*/
    fp_t *aerosol_single_scatter_albedo; /*Aerosol single-scatter albedo (time, column, layer, wavenumber).*/
    fp_t *aerosol_asymmetry_factor; /*Aerosol asymmetry factor (time, column, layer, wavenumber).*/
    fp_t *albedo_grid; /*Wavenumber [cm-1] grid for the surface albedo (wavenumber).*/
    size_t albedo_grid_size; /*Length of albedo_grid array.*/
    Cfc_t *cfc; /*GRTCODE cfc ids and paths.*/
    fp_t **cfc_ppmv; /*CFC abundance [ppmv] (CFC, time, column, level).*/
    fp_t *cloud_fraction; /*Cloud fraction (time, column, layer).*/
    Cia_t *cia; /*GRTCODE CIA ids and paths.*/
    fp_t **cia_ppmv; /*CIA abundance [ppmv] (molecule, time, column, level).*/
    int *cia_species; /*GRTCODE CIA ids.*/
    int clean; /*Flag for clean (aerosol-free) sky.*/
    int clear; /*Flag for clear (cloud-free) sky.*/
    fp_t *emissivity_grid; /*Wavenumber [cm-1] grid for the surface emissivity (wavenumber).*/
    size_t emissivity_grid_size; /*Length of emissivity_grid array.*/
    char h2o_ctm[valuelen]; /*Path to water vapor continuum directory.*/
    fp_t *ice_water_content; /*Cloud ice water content [g m-3] (time, column, layer).*/
    fp_t *layer_pressure; /*Pressure [atm] (time, column, layer).*/
    fp_t *layer_temperature; /*Temperature [K] (time, column, layer).*/
    fp_t *level_pressure; /*Pressure [atm] (time, column, level).*/
    fp_t *level_temperature; /*Temperature [K] (time, column, level).*/
    fp_t *layer_thickness; /*Thickness [m] (time, column, layer).*/
    fp_t *liquid_water_droplet_radius; /*Liquid water equivalent radius [microns] (time, column, layer).*/
    fp_t *liquid_water_content; /*Cloud liquid water content [g m-3] (time, column, layer).*/
    int *molecules; /*GRTCODE molecule ids.*/
    int num_cfcs; /*Length of cfc array.*/
    int num_cia_species; /*Length of cia_species array.*/
    int num_cias; /*Length of cia array.*/
    int num_columns; /*Number of columns.*/
    int num_layers; /*Number of atmospheric layers.*/
    int num_levels; /*Number of atmospheric levels.*/
    int num_molecules; /*Length of molecules array.*/
    int num_times; /*Number of times.*/
    char o3_ctm[valuelen]; /*Path to ozone continuum file.*/
    fp_t **ppmv; /*Molecular abundance [ppmv] (molecule, time, column, level).*/
    fp_t *solar_zenith_angle; /*Cosine of solar zenith angle (time, column).*/
    fp_t *surface_albedo; /*Surface albedo (time, column, wavenumber).*/
    fp_t *surface_emissivity; /*Surface emissivity (time, column, wavenumber).*/
    fp_t *surface_temperature; /*Surface temperature [K] (time, column).*/
    fp_t *total_solar_irradiance; /*Total solar irradiance [W m-2] at TOA (time, column).*/

    int x; /*Starting index for the x-dimension.*/
    int X; /*Ending index for the x-dimension.*/
} Atmosphere_t;


/*The application developer must define this type.*/
typedef struct Output Output_t; /*Structure used to create/write the output.*/


/*The application developer must define this function.*/
void close_flux_file(
    Output_t * const output
); /*Function to close the output dataset.*/


/*The application developer must define this function.*/
void create_flux_file(
    Output_t ** output,
    char const * const path,
    Atmosphere_t const * const atm,
    SpectralGrid_t const * const lw_grid,
    SpectralGrid_t const * const sw_grid,
    int const user_level,
    int const integrated /*Write out integrated flux output.*/
);


Atmosphere_t create_atmosphere(
    Parser_t * const parser
);


void destroy_atmosphere(
    Atmosphere_t * atm
);


void write_output(
    Output_t * output,
    Variables_t id,
    fp_t const * data,
    int time,
    int column
);


#endif
