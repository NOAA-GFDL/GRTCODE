#include <stdlib.h>
#include "cloud_pade_optics.h"
#include "clouds_lib.h"
#include "incomplete_beta.h"
#include "optics_utils.h"
#include "stochastic_clouds.h"


static IncompleteBeta_t beta;
static ty_cloud_optics ice;
static OpticalProperties_t ice_optics;
static ty_cloud_optics liquid;
static OpticalProperties_t liquid_optics;
static TotalWaterPDF_t water_pdf;


/* @brief Set up library global derived types.*/
int initialize_clouds_lib(char const * beta_path, char const * ice_path,
                          char const * liquid_path)
{
    construct_beta(&beta, beta_path);
    construct_cloud_optics(&ice, ice_path);
    construct_optics(&ice_optics, ice.nbnd, ice.band_lims_wvn);
    construct_cloud_optics(&liquid, liquid_path);
    construct_optics(&liquid_optics, liquid.nbnd, liquid.band_lims_wvn);
    construct_water_pdf(&water_pdf, 5, 5, &beta);
    return 0;
}


/* @brief Destroy library global derived types.*/
int finalize_clouds_lib()
{
    destruct_beta(&beta);
    finalize_ty_cloud_optics(&ice);
    destruct_optics(&ice_optics);
    finalize_ty_cloud_optics(&liquid);
    destruct_optics(&liquid_optics);
    destruct_water_pdf(&water_pdf);
    return 0;
}

static double ice_particle_size(double const temperature)
{
    double const tfreeze = 273.16;
    if (temperature > tfreeze - 25.)
    {
        return 100.6;
    }
    else if (temperature > tfreeze - 30.)
    {
        return 80.8;
    }
    else if (temperature > tfreeze - 35.)
    {
        return 93.5;
    }
    else if (temperature > tfreeze - 40.)
    {
        return 63.9;
    }
    else if (temperature > tfreeze - 45.)
    {
        return 42.5;
    }
    else if (temperature > tfreeze - 50.)
    {
        return 39.9;
    }
    else if (temperature > tfreeze - 55.)
    {
        return 21.6;
    }
    else
    {
        return 20.2;
    }
}

int cloud_optics(const double *wavenum,  
                 int const num_wavenum,
                 int const num_layers,
                 const double *mean_cloud_fraction,
                 const double *mean_liquid_content, 
                 const double *mean_ice_content,
                 const double *overlap, 
                 fp_t liquid_radius,
                 const double *temperature, 
                 double *beta_liquid, double *omega_liquid, double *g_liquid,
                 double *beta_ice, double *omega_ice, double *g_ice)
{
    /* Allocate storage for content for each layer */
    double ice_content[num_layers];
    double liquid_content[num_layers];
    /* Presumably, ice_particle_size returns a single value, so we call it once per temperature value. */
    double ice_radius[num_layers];  
    int i, iband;
    for (i = 0; i < num_layers; ++i)
    {
        ice_radius[i] = ice_particle_size(temperature[i])/2.0;
    }

    for (iband = 0; iband < liquid.nbnd; ++iband)
    {
    sample_condensate(water_pdf, num_layers, mean_cloud_fraction, mean_liquid_content,
                      mean_ice_content, overlap, liquid_content, ice_content);
    /* Sample condensate. NOTE: water_pdf is undefined so assume it is defined elsewhere (or remove it) */
        for (i = 0; i < num_layers; ++i)
        {
        int offset = i * num_wavenum;  // offset into the output arrays for this layer
    
            /* Compute optical properties for the liquid case. 
               Note: adjust the arguments to match the actual function signature.
               For example, you might need to provide grid_size and nbnd as well. */

            compute_all_from_pade(&liquid, 
                                  liquid_content[i], liquid_radius, 
                                  &liquid_optics, iband);

            /*Convert from extinction coefficient to optical depth.*/
            map_band_wave(liquid_optics, &iband, wavenum, offset, num_wavenum,
                          beta_liquid, omega_liquid, g_liquid);

            if ((mean_liquid_content[i] > 0) && (mean_cloud_fraction[i] > 0.5) && (iband==1000 || iband==1001) && (i==25))  
            { 
                        fprintf(stdout, "mean_liquid_content = %e\n", mean_liquid_content[i]);  
                        fprintf(stdout, "liquid_optics.extinction_coefficient[%d] = %e\n",
                                iband, liquid_optics.extinction_coefficient[1000]);}                 

            /* Similarly for ice optical properties.
               If ice optical properties depend on ice_rad (calculated above) and ice_content[i]. */
            compute_all_from_pade(&ice, 
                                  ice_content[i], ice_radius[i], 
                                  &ice_optics, iband);
                                  
            map_band_wave(ice_optics, &iband, wavenum, offset, num_wavenum,
                          beta_ice, omega_ice, g_ice);
        }
    }
    return 0;
}

/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
int calculate_overlap(int const num_layers, double const * altitude,
                      double const scale_length, double * alpha)
{
    overlap_parameter(num_layers, altitude, scale_length, alpha);
    return 0;
}