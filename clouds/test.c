#include "cloud_optics.h"
#include "grtcode_utilities.h"


int main(void)
{
    /*Increase verbosity.*/
    grtcode_set_verbosity(4);

    /*Initialize a device.*/
    Device_t device;
    catch(create_device(&device, NULL));

    /*Set up a test spectral grid.*/
    SpectralGrid_t grid;
    catch(create_spectral_grid(&grid, 1., 3250., 0.1));

    /*Make up some test data.*/
    int const num_layers = 5;
    fp_t pressure[num_layers] = {1., 10., 100., 500., 1000.}; /*[mb].*/
    fp_t temperature[num_layers] = {230., 235., 250., 275., 295,}; /*[K].*/
    fp_t cloud_fraction[num_layers] = {0.65, 0.23, 0.11, 0.75, 0.8};
    fp_t liquid_content[num_layers] = {1.e-6, 5.e-7, 3.e-6, 2.e-7, 1.1e-7};
    fp_t ice_content[num_layers] = {5.12e-7, 4.4e-7, 1.2e-6, 9.2e-7, 2.1e-6};
    fp_t thickness[num_layers] = {3., 1.5, 6., 9., 7.77};

    char * beta_path = "../grtcode-data/clouds/beta_distribution.nc";
    char * liquid_path = "../grtcode-data/clouds/hu_stamnes.nc";
    char * ice_path = "../grtcode-data/clouds/fu_liou.nc";

    /*Initialize the cloud optics.*/
    CloudOptics_t cloud_optics;
    catch(create_cloud_optics(&cloud_optics, beta_path, liquid_path, ice_path,
                              num_layers, grid, device));

    /*Calculate the cloud optics.*/
    int j;
    for (j=0; j<10000; ++j)
    {
        printf("start: iteration %d\n", j);
        catch(calculate_cloud_optics(&cloud_optics, num_layers, pressure, temperature,
                                     thickness, cloud_fraction, liquid_content, ice_content));
    }

    /*Clean up.*/
    catch(destroy_cloud_optics(&cloud_optics));
    return 0;
}
