/* @brief Ice water cloud optics parameterizations from doi: 10.1175/2009JCLI2844.1*/
#include <stdlib.h>
#include "debug.h"
#include "ice_cloud_optics.h"
#include "netcdf.h"
#include "netcdf_utils.h"
#include "optics_utils.h"


/* @brief Donner parameterization for ice particle size.*/
HOST DEVICE static fp_t ice_particle_size(
    fp_t const temperature /**< Temperature [K].*/
)
{
    fp_t const tfreeze = 273.16;
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


/* @brief Calculates cloud optics.*/
HOST DEVICE static void optics(
    int const num_radius_bins, /**< Number of radius bins.*/
    fp_t const * radii, /**< Radii [micron] (radius).*/
    int const num_order, /**< Number of polynomial orders.*/
    int const num_bands, /**< Number of bands.*/
    int const last_ir_band, /**< Index of the last infrared band.*/
    fp_t const * a, /**< a parameter (band, order).*/
    fp_t const * b, /**< b parameter (band, order).*/
    fp_t const * c, /**< c parameter (band, order).*/
    fp_t const ice_concentration, /**< Ice concentration [g m-3].*/
    fp_t const equivalent_radius, /**< Equivalent radius [micron].*/
    fp_t const scale_factor, /**< Scaling factor.*/
    fp_t const temperature, /**< Temperature [K].*/
    fp_t const thickness, /**< Layer thickness [m].*/
    int const band, /**< Band index.*/
    fp_t * optical_depth, /**< Optical depth.*/
    fp_t * single_scatter_albedo, /**< Single-scatter albedo.*/
    fp_t * asymmetry_factor /**< Asymmetry factor.*/
)
{
    if (ice_concentration > 0.)
    {
        fp_t radius = equivalent_radius > 0. ? equivalent_radius : ice_particle_size(temperature);
        fp_t const min_radius = 13.;
        radius = min_radius > scale_factor*radius ? min_radius : scale_factor*radius;
        int r;
        for (r=1; r<num_radius_bins; ++r)
        {
            if (radii[2*r] > radius) break;
        }
        r -= 1;
        fp_t d[16];
        d[0] = 1.;
        int i;
        for (i=1; i<num_order; ++i)
        {
            d[i] = d[i - 1]*radius;
        }
        fp_t d_inv[16];
        for (i=0; i<num_order; ++i)
        {
            d_inv[i] = 1./d[i];
        }

        fp_t asum = 0.;
        for (i=0; i<num_order; ++i)
        {
            asum += a[band*num_order + i]*d_inv[i];
        }
        fp_t const extinction = ice_concentration*asum;
        *optical_depth = extinction*thickness;
        if (band < last_ir_band)
        {
            fp_t bsum = 0.;
            for (i=0; i<num_order; ++i)
            {
                bsum += b[band*num_order + i]*d_inv[i];
            }
            *single_scatter_albedo = 1. - (ice_concentration*bsum/(extinction));
        }
        else
        {
            fp_t bsum = 0.;
            for (i=0; i<num_order; ++i)
            {
                bsum += b[band*num_order + i]*d[i];
            }
            *single_scatter_albedo = 1. - bsum;
        }
        fp_t csum = 0.;
        for (i=0; i<num_order; ++i)
        {
            csum += c[(r*num_bands + band)*num_order + i]*d[i];
        }
        *asymmetry_factor = csum;
    }
    else
    {
        *optical_depth = 0.;
        *single_scatter_albedo = 0.;
        *asymmetry_factor = 0.;
    }
    return;
}


/* @brief Constructs a IceCloudOptics object.*/
EXTERN int create_ice_optics(IceCloudOptics_t * self, char const * path,
                             Device_t const * const device)
{
    not_null(self);
    not_null(path);
    not_null(device);
    self->device = *device;
    int ncid;
    catch(open_dataset(path, &ncid));
    catch(read_dimlen(ncid, "band", &(self->num_bands)));
    catch(read_dimlen(ncid, "order5", &(self->num_order)));
    catch(read_dimlen(ncid, "radius", &(self->num_radius_bins)));
    nc_type var_type;
    if (sizeof(fp_t) == sizeof(double))
    {
        var_type = NC_DOUBLE;
    }
    else
    {
        var_type = NC_FLOAT;
    }
    catch(read_variable(ncid, "radius_bnds", (void **)&(self->radii), var_type,
                        self->device, NULL, NULL));
    fp_t * band_limits;
    catch(read_variable(ncid, "band_bnds", (void **)(&band_limits), var_type, HOST_ONLY,
                        NULL, NULL));
    if (self->device == HOST_ONLY)
    {
        self->band_limits = band_limits;
    }
    else
    {
        gmalloc(self->band_limits, 2*self->num_bands, self->device);
        gmemcpy(self->band_limits, band_limits, 2*self->num_bands,
                self->device, FROM_HOST);
    }
    fp_t * buffer;
    buffer = (fp_t *)malloc(sizeof(fp_t)*self->num_bands);
    int i;
    for (i=0; i<self->num_bands; ++i)
    {
        buffer[i] = 0.5*(band_limits[2*i] + band_limits[2*i + 1]);
    }
    if (self->device == HOST_ONLY)
    {
        self->bands = buffer;
    }
    else
    {
        gmalloc(self->bands, self->num_bands, self->device);
        gmemcpy(self->bands, buffer, self->num_bands, self->device, FROM_HOST);
        free(band_limits);
        free(buffer);
    }
    catch(read_attribute(ncid, "band_bnds", "last_IR_band", NC_INT, &(self->last_ir_band)));
    catch(read_variable(ncid, "a", (void **)&(self->a), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "b", (void **)&(self->b), var_type, self->device, NULL, NULL));
    catch(read_variable(ncid, "c", (void **)&(self->c), var_type, self->device, NULL, NULL));
    catch(close_dataset(ncid));
    return GRTCODE_SUCCESS;
}


/* @brief Constructs a IceCloudOptics object.*/
EXTERN int destroy_ice_optics(IceCloudOptics_t * self)
{
    gfree(self->radii, self->device);
    gfree(self->band_limits, self->device);
    gfree(self->bands, self->device);
    gfree(self->a, self->device);
    gfree(self->b, self->device);
    gfree(self->c, self->device);
    return GRTCODE_SUCCESS;
}


/* @brief Calculates cloud optics.*/
int calculate_ice_optics(int const num_radius_bins, fp_t const * radii,
                         int const num_order, int const num_bands, int const last_ir_band,
                         fp_t const * a, fp_t const * b, fp_t const * c,
                         int const num_layers,
                         fp_t const * ice_concentration,
                         fp_t const equivalent_radius,
                         fp_t const scale_factor,
                         fp_t const * temperature,
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
        optics(num_radius_bins, radii, num_order, num_bands, last_ir_band, a, b, c,
               ice_concentration[i], equivalent_radius, scale_factor,
               temperature[layer], thickness[layer], band, &(optical_depth[i]),
               &(single_scatter_albedo[i]), &(asymmetry_factor[i]));
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/* @brief Calculates cloud optics.*/
__global__ void calculate_ice_optics_d(int const num_radius_bins, fp_t const * radii,
                                       int const num_order, int const num_bands,
                                       int const last_ir_band,
                                       fp_t const * a, fp_t const * b, fp_t const * c,
                                       int const num_layers,
                                       fp_t const * ice_concentration,
                                       fp_t const equivalent_radius,
                                       fp_t const scale_factor,
                                       fp_t const * temperature,
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
        optics(num_radius_bins, radii, num_order, num_bands, last_ir_band, a, b, c,
               ice_concentration[i], equivalent_radius, scale_factor,
               temperature[layer], thickness[layer], band, &(optical_depth[i]),
               &(single_scatter_albedo[i]), &(asymmetry_factor[i]));
    }
    return;
}
#endif
