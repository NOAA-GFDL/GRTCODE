#include <math.h>
#include "cloud_optics.h"
#include "cloud_utils.h"
#include "debug.h"
#include "device.h"
#include "ice_cloud_optics.h"
#include "incomplete_beta.h"
#include "liquid_cloud_optics.h"
#include "optics.h"
#include "spectral_grid.h"
#include "stochastic_clouds.h"


int create_cloud_optics(CloudOptics_t * self, char const * beta_path, char const * liquid_path,
                        char const * ice_path, int const num_layers, SpectralGrid_t const grid,
                        Device_t const device)
{
    not_null(self);
    not_null(beta_path);
    not_null(liquid_path);
    not_null(ice_path);
    self->device = device;
    self->num_layers = num_layers;
    catch(create_spectral_grid(&(self->grid), grid.w0, grid.wn, grid.dw));

    /*Copy cloud parameterization input data to the device.*/
    catch(create_beta(&(self->beta), beta_path, &device));
    catch(create_water_pdf(&(self->pdf), 5, 5, &(self->beta)));
    catch(create_liquid_optics(&(self->liquid), liquid_path, &device));
    catch(create_ice_optics(&(self->ice), ice_path, &device));
    self->num_subcolumns = self->liquid.num_bands > self->ice.num_bands ?
                           self->liquid.num_bands : self->ice.num_bands;

    /*Create mapping between bands and grid.*/
    gmalloc(self->liquid_mapping, grid.n, device);
    glaunch(grid_band_mapping, grid.n, device, grid.n, grid.w0, grid.dw, self->liquid.num_bands,
            self->liquid.band_limits, self->liquid_mapping);
    gmalloc(self->ice_mapping, grid.n, device);
    glaunch(grid_band_mapping, grid.n, device, grid.n, grid.w0, grid.dw, self->ice.num_bands,
            self->ice.band_limits, self->ice_mapping);

    /*Allocate optics arrays.*/
    gmalloc(self->tau, num_layers*self->num_subcolumns, device);
    gmalloc(self->omega, num_layers*self->num_subcolumns, device);
    gmalloc(self->g, num_layers*self->num_subcolumns, device);
    self->random = (fp_t *)malloc(sizeof(fp_t)*num_layers*self->num_subcolumns);
    catch(create_optics(&(self->liquid_optics), num_layers, &grid, &device));
    catch(create_optics(&(self->ice_optics), num_layers, &grid, &device));

    /*Allocate cloud properties arrays.*/
    gmalloc(self->temperature, num_layers, device);
    gmalloc(self->cc, num_layers, device);
    gmalloc(self->clwc, num_layers, device);
    gmalloc(self->ciwc, num_layers, device);
    gmalloc(self->height, num_layers, device);
    gmalloc(self->thickness, num_layers, device);
    gmalloc(self->alpha, num_layers-1, device);
    gmalloc(self->random1_d, num_layers*self->num_subcolumns, device);
    gmalloc(self->random2_d, (num_layers-1)*self->num_subcolumns, device);
    gmalloc(self->liquid_condensate, num_layers*self->num_subcolumns, device);
    gmalloc(self->ice_condensate, num_layers*self->num_subcolumns, device);
    return GRTCODE_SUCCESS;
}


int destroy_cloud_optics(CloudOptics_t * self)
{
    catch(destroy_ice_optics(&(self->ice)));
    catch(destroy_liquid_optics(&(self->liquid)));
    catch(destroy_water_pdf(&(self->pdf)));
    catch(destroy_beta(&(self->beta)));
    catch(destroy_optics(&(self->liquid_optics)));
    catch(destroy_optics(&(self->ice_optics)));
    gfree(self->tau, self->device);
    gfree(self->omega, self->device);
    gfree(self->g, self->device);
    gfree(self->temperature, self->device);
    gfree(self->cc, self->device);
    gfree(self->clwc, self->device);
    gfree(self->ciwc, self->device);
    gfree(self->height, self->device);
    gfree(self->thickness, self->device);
    gfree(self->alpha, self->device);
    gfree(self->liquid_condensate, self->device);
    gfree(self->ice_condensate, self->device);
    gfree(self->random1_d, self->device);
    gfree(self->random2_d, self->device);
    gfree(self->liquid_mapping, self->device);
    gfree(self->ice_mapping, self->device);
    free(self->random);
    return GRTCODE_SUCCESS;
}


int calculate_cloud_optics(CloudOptics_t * self, int const num_layers,
                           fp_t const * pressure, fp_t const * temperature,
                           fp_t const * thickness, fp_t const * cloud_fraction,
                           fp_t const * liquid_content, fp_t const * ice_content)
{
    /*Defense.*/
    not_null(self);
    not_null(pressure);
    not_null(temperature);
    not_null(thickness);
    not_null(cloud_fraction);
    not_null(liquid_content);
    not_null(ice_content);
    if (num_layers > self->num_layers)
    {
        raise(GRTCODE_VALUE_ERR, "%s", "not enough space allocated on the device.");
    }

    /*Calculate the altitude.*/
    fp_t const equivalent_radius = 10.; /*[microns].*/
    fp_t const scale_length = 2.; /*[km].*/
    fp_t const mb_to_Pa = 100.;
    fp_t altitude[num_layers];
    int i;
    for (i=0; i<num_layers; ++i)
    {
        altitude[i] = log(mb_to_Pa*pressure[i])*7.3; /*[km].*/
    }

    /*Copy cloud properties to the device.*/
    gmemcpy(self->temperature, temperature, num_layers, self->device, FROM_HOST);
    gmemcpy(self->cc, cloud_fraction, num_layers, self->device, FROM_HOST);
    gmemcpy(self->clwc, liquid_content, num_layers, self->device, FROM_HOST);
    gmemcpy(self->ciwc, ice_content, num_layers, self->device, FROM_HOST);
    gmemcpy(self->height, altitude, num_layers, self->device, FROM_HOST);
    gmemcpy(self->thickness, thickness, num_layers, self->device, FROM_HOST);

    /*Calculate overlap parameters.*/
    glaunch(overlap_parameter, num_layers-1, self->device, num_layers, self->height,
            scale_length, self->alpha);

    /*For each layer of each subcolumn, generate some random numbers.*/
    /*Random numbers for total water vapor sampling.*/
    for (i=0; i<num_layers; ++i)
    {
        self->random[i] = ((fp_t)rand())/((fp_t)RAND_MAX);
    }
    /*Copy the first subcolumn for all the rest of the subcolumns.*/
    int j;
    for (j=1; j<self->num_subcolumns; ++j)
    {
        for (i=0; i<num_layers; ++i)
        {
            self->random[j*num_layers + i] = self->random[i];
        }
    }
    gmemcpy(self->random1_d, self->random, num_layers*self->num_subcolumns,
            self->device, FROM_HOST);

    /*Random numbers for overlap (correlation between layers).*/
    for (i=0; i<(num_layers-1); ++i)
    {
        self->random[i] = ((fp_t)rand())/((fp_t)RAND_MAX);
    }
    /*Copy the first subcolumn for all the rest of the subcolumns.*/
    for (j=1; j<self->num_subcolumns; ++j)
    {
        for (i=0; i<(num_layers-1); ++i)
        {
            self->random[j*(num_layers-1) + i] = self->random[i];
        }
    }
    gmemcpy(self->random2_d, self->random, (num_layers-1)*self->num_subcolumns,
            self->device, FROM_HOST);

    /*For each layer, sample the condensate once per each subcolumn*/
    glaunch(sample_condensate, self->num_subcolumns, self->device,
            self->pdf.num_shape, self->pdf.num_x, self->pdf.x, self->pdf.y,
            self->pdf.y_inverse, self->pdf.p, self->pdf.q,
            self->num_subcolumns, num_layers, self->cc, self->clwc, self->ciwc,
            self->alpha, self->random1_d, self->random2_d, self->liquid_condensate,
            self->ice_condensate);

    /*For each layer, calculate the liquid cloud optics for each subcolumn.*/
    int num_bands = self->liquid.num_bands;
    glaunch(calculate_liquid_optics, num_layers*num_bands, self->device,
            self->liquid.min_radius, self->liquid.max_radius, self->liquid.num_radius_bins,
            self->liquid.radii, num_bands, self->liquid.a1, self->liquid.b1,
            self->liquid.c1, self->liquid.a2, self->liquid.b2, self->liquid.c2,
            self->liquid.a3, self->liquid.b3, self->liquid.c3, num_layers,
            self->liquid_condensate, equivalent_radius, self->thickness,
            self->tau, self->omega, self->g);

    /*Expand the optics from bands to spectral grid.*/
    glaunch(process_optics, self->grid.n, self->device, self->grid.n, num_layers, num_bands,
            self->liquid_mapping, self->tau, self->omega, self->g, self->liquid_optics.tau,
            self->liquid_optics.omega, self->liquid_optics.g);

    /*For each layer, calculate the ice cloud optics for each subcolumn.*/
    num_bands = self->ice.num_bands;
    fp_t const radius = -1.;
    fp_t const scale_factor = 1.;
    glaunch(calculate_ice_optics, num_layers*num_bands, self->device,
            self->ice.num_radius_bins, self->ice.radii, self->ice.num_order, num_bands,
            self->ice.last_ir_band, self->ice.a, self->ice.b, self->ice.c, num_layers,
            self->ice_condensate, radius, scale_factor, self->temperature, self->thickness,
            self->tau, self->omega, self->g);

    /*Expand the optics from bands to spectral grid.*/
    glaunch(process_optics, self->grid.n, self->device, self->grid.n, num_layers, num_bands,
            self->ice_mapping, self->tau, self->omega, self->g, self->ice_optics.tau,
            self->ice_optics.omega, self->ice_optics.g);
    return GRTCODE_SUCCESS;
}
