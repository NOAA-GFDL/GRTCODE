#include "cloud_utils.h"
#include "debug.h"


/** @brief Create a map from wavenumber to cloud band.*/
int grid_band_mapping(int const grid_size, fp_t const w0,
                      fp_t const dw, int const num_bands,
                      fp_t * const band_limits, int * mapping)
{
    int i;
    for (i=0; i<grid_size; ++i)
    {
        int band;
        fp_t const w = w0 + i*dw;
        if (w < band_limits[0])
        {
            band = 0;
        }
        else if (w > band_limits[2*num_bands-1])
        {
            band = num_bands - 1;
        }
        else
        {
            for (band=0; band<num_bands-1; ++band)
            {
                if (w >= band_limits[band*2] && w <= band_limits[band*2+1])
                {
                    break;
                }
            }
        }
        mapping[i] = band;
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/** @brief Create a map from wavenumber to cloud band.*/
__global__ void grid_band_mapping_d(int const grid_size, fp_t const w0,
                                    fp_t const dw, int const num_bands,
                                    fp_t * const band_limits, int * mapping)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < grid_size)
    {
        int band;
        fp_t const w = w0 + i*dw;
        if (w < band_limits[0])
        {
            band = 0;
        }
        else if (w > band_limits[2*num_bands-1])
        {
            band = num_bands - 1;
        }
        else
        {
            for (band=0; band<num_bands-1; ++band)
            {
                if (w >= band_limits[band*2] && w <= band_limits[band*2+1])
                {
                    break;
                }
            }
        }
        mapping[i] = band;
    }
    return;
}
#endif


int process_optics(int const grid_size, int const num_layers, int const num_bands,
                   int const * mapping, fp_t const * tau, fp_t const * omega,
                   fp_t const * g, fp_t * optics_tau,
                   fp_t * optics_omega, fp_t * optics_g)
{
    int j;
    for (j=0; j<num_layers; ++j)
    {
        int i;
        for (i=0; i<grid_size; ++i)
        {
            optics_tau[j*grid_size + i] = tau[mapping[i]*num_layers + j];
            optics_omega[j*grid_size + i] = omega[mapping[i]*num_layers + j];
            optics_g[j*grid_size + i] = g[mapping[i]*num_layers + j];
        }
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
__global__ void process_optics_d(int const grid_size, int const num_layers, int const num_bands,
                                 int const * mapping, fp_t const * tau, fp_t const * omega,
                                 fp_t const * g, fp_t * optics_tau,
                                 fp_t * optics_omega, fp_t * optics_g)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < grid_size)
    {
        int j;
        for (j=0; j<num_layers; ++j)
        {
            optics_tau[j*grid_size + i] = tau[mapping[i]*num_layers + j];
            optics_omega[j*grid_size + i] = omega[mapping[i]*num_layers + j];
            optics_g[j*grid_size + i] = g[mapping[i]*num_layers + j];
        }
    }
    return;
}
#endif
