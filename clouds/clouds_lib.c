#include <stdlib.h>
#include "hu_stamnes.h"
#include "ice_cloud_optics.h"
#include "incomplete_beta.h"
#include "optics_utils.h"
#include "stochastic_clouds.h"


static IncompleteBeta_t beta;
static IceCloudOptics_t ice;
static OpticalProperties_t ice_optics;
static HuStamnes_t liquid;
static OpticalProperties_t liquid_optics;
static TotalWaterPDF_t water_pdf;


/* @brief Set up library global derived types.*/
int initialize_clouds_lib(char const * beta_path, char const * ice_path,
                          char const * liquid_path, int const * beta_shape)
{
    construct_beta(&beta, beta_path);
    construct_ice_optics(&ice, ice_path);
    construct_optics(&ice_optics, ice.num_bands, ice.bands, ice.band_limits);
    construct_liquid_optics(&liquid, liquid_path);
    construct_optics(&liquid_optics, liquid.num_bands, liquid.bands, liquid.band_limits);
    int beta_shape_ = beta_shape == NULL ? 5 : *beta_shape;
    construct_water_pdf(&water_pdf, beta_shape_, beta_shape_, &beta);
    return 0;
}


/* @brief Destroy library global derived types.*/
int finalize_clouds_lib()
{
    destruct_beta(&beta);
    destruct_ice_optics(&ice);
    destruct_optics(&ice_optics);
    destruct_liquid_optics(&liquid);
    destruct_optics(&liquid_optics);
    destruct_water_pdf(&water_pdf);
    return 0;
}


/* @brief Sample ice and liquid cloud optical properties for a column.*/
int cloud_optics(int const num_bands, double const * band_centers, double const * band_limits,
                 int const num_layers, double const * mean_cloud_fraction,
                 double const * mean_liquid_content, double const *mean_ice_content,
                 double const * overlap, double const liquid_radius,
                 double const * temperature, double * beta_liquid, double * omega_liquid,
                 double * g_liquid, double * beta_ice, double *omega_ice, double * g_ice)
{
    OpticalProperties_t remapped_ice_optics;
    construct_optics(&remapped_ice_optics, num_bands, band_centers, band_limits);
    OpticalProperties_t remapped_liquid_optics;
    construct_optics(&remapped_liquid_optics, num_bands, band_centers, band_limits);

    double ice_content[num_layers];
    double liquid_content[num_layers];
    sample_condensate(water_pdf, num_layers, mean_cloud_fraction, mean_liquid_content,
                      mean_ice_content, overlap, liquid_content, ice_content);

    int i;
    for (i=0; i<num_layers; ++i)
    {
        int offset = i*num_bands;
        if (liquid_content[i] > 0.)
        {
            calculate_liquid_optics(liquid, liquid_content[i], liquid_radius,
                                    &liquid_optics);
            thick_average(liquid_optics, &remapped_liquid_optics, NULL,
                          &liquid.last_ir_band);
            int j;
            for (j=0; j<num_bands; ++j)
            {
                beta_liquid[offset + j] = remapped_liquid_optics.extinction_coefficient[j];
                omega_liquid[offset + j] = remapped_liquid_optics.single_scatter_albedo[j];
                g_liquid[offset + j] = remapped_liquid_optics.asymmetry_factor[j];
            }
        }
        else
        {
            int j;
            for (j=0; j<num_bands; ++j)
            {
                beta_liquid[offset + j] = 0.;
                omega_liquid[offset + j] = 0.;
                g_liquid[offset + j] = 0.;
            }
        }
        if (ice_content[i] > 0.)
        {
            calculate_ice_optics(ice, ice_content[i], -1., 1., temperature[i], &ice_optics);
            thick_average(ice_optics, &remapped_ice_optics, NULL, &ice.last_ir_band);
            int j;
            for (j=0; j<num_bands; ++j)
            {
                beta_ice[offset + j] = remapped_ice_optics.extinction_coefficient[j];
                omega_ice[offset + j] = remapped_ice_optics.single_scatter_albedo[j];
                g_ice[offset + j] = remapped_ice_optics.asymmetry_factor[j];
            }
        }
        else
        {
            int j;
            for (j=0; j<num_bands; ++j)
            {
                beta_ice[offset + j] = 0.;
                omega_ice[offset + j] = 0.;
                g_ice[offset + j] = 0.;
            }
        }
    }
    destruct_optics(&remapped_liquid_optics);
    destruct_optics(&remapped_ice_optics);
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
