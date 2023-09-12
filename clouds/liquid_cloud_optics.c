/* @brief Hu and Stamnes cloud liquid water parameterization from
          doi: 10.1175/1520-0442(1993)006<0728:AAPOTR>2.0.CO;2*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "liquid_cloud_optics.h"
#include "netcdf.h"
#include "netcdf_utils.h"
#include "optics_utils.h"


/* @brief Calculates cloud optics for a single band.*/
HOST DEVICE static void optics(
    fp_t const min_radius, /**< Minimum radius [micron].*/
    fp_t const max_radius, /**< Maximum radius [micron].*/
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_bands, /**< Number of bands.*/
    fp_t const * a1, /**< a1 parameter (radius, band).*/
    fp_t const * b1, /**< b1 parameter (radisu, band).*/
    fp_t const * c1, /**< c1 parameter (radius, band).*/
    fp_t const * a2, /**< a2 parameter (radius, band).*/
    fp_t const * b2, /**< b2 parameter (radius, band).*/
    fp_t const * c2, /**< c2 parameter (radius, band).*/
    fp_t const * a3, /**< a3 parameter (radius, band).*/
    fp_t const * b3, /**< b3 parameter (radius, band).*/
    fp_t const * c3, /**< c3 parameter (radius, band).*/
    fp_t const water_concentration, /**< Water concentration [g m-3].*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const thickness, /**< Layer thickness [m].*/
    int const band, /**< Band index.*/
    fp_t * optical_depth, /**< Optical depth.*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo.*/
    fp_t * asymmetry_factor /**< Asymmetry factor.*/
)
{
    if (water_concentration > 0.)
    {
        fp_t r = min_radius > equivalent_radius ? min_radius : equivalent_radius;
        r = max_radius < r ? max_radius : r;
        int i;
        for (i=1; i<num_radius_bins; ++i)
        {
            if (radii[i] > r) break;
        }
        i = (i - 1)*num_bands + band;
        fp_t const m_to_km = 1.e-3; /*[km m-1].*/
        *optical_depth = water_concentration*m_to_km*(a1[i]*
                         pow(r, b1[i]) + c1[i])*thickness; /*Equation 13.*/
        *single_scatter_albedo = 1. - (a2[i]*pow(r, b2[i]) + c2[i]); /*Equation 14.*/
        if (*single_scatter_albedo < 0.)
        {
            *single_scatter_albedo = 0.;
        }
        *asymmetry_factor = a3[i]*pow(r, b3[i]) + c3[i]; /*Equation 15.*/
    }
    else
    {
        *optical_depth = 0.;
        *single_scatter_albedo = 0.;
        *asymmetry_factor = 0.;
    }
    return;
}


/* @brief Constructs a LiquidCloudOptics object.*/
EXTERN int create_liquid_optics(LiquidCloudOptics_t * self, char const * path,
                                Device_t const * const device)
{
    not_null(self);
    not_null(path);
    not_null(device);
    self->device = *device;
    int ncid;
    catch(open_dataset(path, &ncid));
    char name[128];
    memset(name, '\0', 128);
    catch(read_attribute(ncid, "radius", "bounds", NC_CHAR, name));
    float radii_range[2];
    catch(read_attribute(ncid, "radius", "valid_range", NC_FLOAT, radii_range));
    self->min_radius = (fp_t)(radii_range[0]);
    self->max_radius = (fp_t)(radii_range[1]);
    fp_t * radius_bounds;
    nc_type var_type;
    if (sizeof(fp_t) == sizeof(double))
    {
        var_type = NC_DOUBLE;
    }
    else
    {
        var_type = NC_FLOAT;
    }
    catch(read_variable(ncid, name, (void **)(&radius_bounds), var_type, HOST_ONLY, NULL, NULL));
    catch(read_dimlen(ncid, "band", &(self->num_bands)));
    catch(read_dimlen(ncid, "radius", &(self->num_radius_bins)));
    fp_t * radii = (fp_t *)malloc(sizeof(fp_t)*(self->num_radius_bins + 1));
    int i;
    for (i=0; i<self->num_radius_bins; ++i)
    {
        radii[i] = radius_bounds[2*i];
    }
    radii[self->num_radius_bins] = radius_bounds[2*(self->num_radius_bins - 1) + 1];
    free(radius_bounds);
    if (self->device == HOST_ONLY)
    {
        self->radii = radii;
    }
    else
    {
        gmalloc(self->radii, self->num_radius_bins + 1, self->device);
        gmemcpy(self->radii, radii, self->num_radius_bins + 1, self->device, FROM_HOST);
        free(radii);
    }
    catch(read_variable(ncid, "band_bnds", (void **)&(self->band_limits), var_type,
                        self->device, NULL, NULL));
    catch(read_attribute(ncid, "band_bnds", "last_IR_band", NC_INT, &self->last_ir_band));
    catch(read_variable(ncid, "band", (void **)&(self->bands), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "a1", (void **)&(self->a1), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "a2", (void **)&(self->a2), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "a3", (void **)&(self->a3), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "b1", (void **)&(self->b1), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "b2", (void **)&(self->b2), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "b3", (void **)&(self->b3), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "c1", (void **)&(self->c1), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "c2", (void **)&(self->c2), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "c3", (void **)&(self->c3), var_type, self->device, NULL, NULL));
    catch(close_dataset(ncid));
    return GRTCODE_SUCCESS;
}


/* @brief Destructs a HuStamnes object.*/
EXTERN int destroy_liquid_optics(LiquidCloudOptics_t * self)
{
    gfree(self->a1, self->device);
    gfree(self->a2, self->device);
    gfree(self->a3, self->device);
    gfree(self->band_limits, self->device);
    gfree(self->bands, self->device);
    gfree(self->b1, self->device);
    gfree(self->b2, self->device);
    gfree(self->b3, self->device);
    gfree(self->c1, self->device);
    gfree(self->c2, self->device);
    gfree(self->c3, self->device);
    gfree(self->radii, self->device);
    return GRTCODE_SUCCESS;
}


/* @brief Calculates liquid cloud optics for all bands.*/
int calculate_liquid_optics(fp_t const min_radius, fp_t const max_radius,
                            int const num_radius_bins, fp_t const * radii,
                            int const num_bands, fp_t const * a1, fp_t const * b1,
                            fp_t const * c1, fp_t const * a2, fp_t const * b2,
                            fp_t const * c2, fp_t const * a3, fp_t const * b3,
                            fp_t const * c3, int const num_layers,
                            fp_t const * water_concentration,
                            fp_t const equivalent_radius,
                            fp_t const * thickness,
                            fp_t * optical_depth,
                            fp_t * single_scatter_albedo,
                            fp_t * asymmetry_factor)
{
    int i;
    for (i=0; i<num_bands*num_layers; ++i)
    {
        int const band = i/num_layers;
        int const layer = i - band*num_layers;
        optics(min_radius, max_radius, num_radius_bins, radii, num_bands, a1, b1, c1,
               a2, b2, c2, a3, b3, c3, water_concentration[i],
               equivalent_radius, thickness[layer], band, &(optical_depth[i]),
               &(single_scatter_albedo[i]), &(asymmetry_factor[i]));
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/* @brief Calculates liquid cloud optics for all bands.*/
__global__ void calculate_liquid_optics_d(fp_t const min_radius, fp_t const max_radius,
                                          int const num_radius_bins, fp_t const * radii,
                                          int const num_bands, fp_t const * a1, fp_t const * b1,
                                          fp_t const * c1, fp_t const * a2, fp_t const * b2,
                                          fp_t const * c2, fp_t const * a3, fp_t const * b3,
                                          fp_t const * c3, int const num_layers,
                                          fp_t const * water_concentration,
                                          fp_t const equivalent_radius,
                                          fp_t const * thickness,
                                          fp_t * optical_depth,
                                          fp_t * single_scatter_albedo,
                                          fp_t * asymmetry_factor)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_bands*num_layers)
    {
        int const band = i/num_layers;
        int const layer = i - band*num_layers;
        optics(min_radius, max_radius, num_radius_bins, radii, num_bands, a1, b1, c1,
               a2, b2, c2, a3, b3, c3, water_concentration[i],
               equivalent_radius, thickness[layer], band, &(optical_depth[i]),
               &(single_scatter_albedo[i]), &(asymmetry_factor[i]));
    }
    return;
}
#endif
