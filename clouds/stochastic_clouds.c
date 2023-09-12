/* @brief Subcolumn generator for calculating stochastic clouds.*/
#include <math.h>
#include <stdlib.h>
#include "debug.h"
#include "incomplete_beta.h"
#include "stochastic_clouds.h"


/** @brief Calculates random array for cloudiness from equations 8 and 9 from
           doi: 10.1256/qj.03.99.*/
HOST DEVICE static void cloudiness(
    int const num_layers, /**< Number of layers.*/
    fp_t const * overlap_parameter, /**< Overlap parameter (layer - 1).*/
    fp_t const * random1, /**< Array of random numbers (layer).*/
    fp_t const * random2, /**< Array of random numbers (layer - 1).*/
    fp_t * x /**< Parameter x (layer).*/
)
{
    int i;
    for (i=0; i<num_layers; ++i)
    {
        x[i] = random1[i];
    }
    for (i=0; i<num_layers - 1; ++i)
    {
        if (random2[i] <= overlap_parameter[i])
        {
            x[i+1] = x[i];
        }
    }
    return;
}


/** @brief Calculates normalized saturation specific humidity from equation A1 from
           doi: 10.1175/MWR3257.1*/
HOST DEVICE static fp_t specific_saturation_humidity(
    int const num_shape, /**< Number of beta distribution shape parameters.*/
    int const num_x, /**< Number of beta distribution x values.*/
    fp_t const * x, /**< Beta distribution x values (x).*/
    fp_t const * y_inverse, /**< Beta distribution deviate values (shape, shape, x).*/
    int const p, /**< Beta distribution shape parameter.*/
    int const q, /**< Beta distribution shape parameter.*/
    fp_t const cloud_fraction /**< Grid-cell mean cloud fraction.*/
)
{
    return beta_inverse(num_shape, num_x, x, y_inverse, p, q, 1. - cloud_fraction);
}


/** @brief Calculates the width of the total water probability distribution function (b - a)
           from equation A2 from doi: 10.1175/MWR3257.1, ignoring the parameter alpha.*/
HOST DEVICE static fp_t width(
    int const num_shape, /**< Number of beta distribution shape parameters.*/
    int const num_x, /**< Number of beta distribution x values.*/
    fp_t const * x, /**< Beta distribution x values (x).*/
    fp_t const * y, /**< Beta distribution values (shape, shape, x).*/
    int const p, /**< Beta distribution shape parameter.*/
    int const q, /**< Beta distribution shape parameter.*/
    fp_t const cloud_fraction, /**< Grid-cell mean cloud fraction.*/
    fp_t const lwc, /**< Grid-cell mean liquid water content.*/
    fp_t const iwc, /**< Grid-cell mean ice cloud content.*/
    fp_t const qs /**< Saturation specific humidity.*/
)
{
    fp_t beta = beta_value(num_shape, num_x, x, y, p + 1, q, qs);
    return (lwc + iwc)/((((fp_t)p)/((fp_t)(p + q)))*(1. - beta) - qs*cloud_fraction);
}


/* @brief Constructs a TotalWaterPDF object.*/
EXTERN int create_water_pdf(TotalWaterPDF_t * self, int const p, int const q,
                            IncompleteBeta_t const * beta)
{
    self->device = beta->device;
    self->num_shape = beta->num_shape;
    self->num_x = beta->num_x;
    self->x = beta->x;
    self->y = beta->y;
    self->y_inverse = beta->y_inverse;
    self->p = p;
    self->q = q;
    return GRTCODE_SUCCESS;
}


/* @brief Destructs a TotalWaterPDF object.*/
EXTERN int destroy_water_pdf(TotalWaterPDF_t * self)
{
    self->x = NULL;
    self->y = NULL;
    self->y_inverse = NULL;
    return GRTCODE_SUCCESS;
}


/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
int overlap_parameter(int const num_layers, fp_t const * altitude,
                      fp_t const scale_length, fp_t * alpha)
{
    int i;
    for (i=0; i<num_layers - 1; ++i)
    {
        alpha[i] = exp(-1.*fabs(altitude[i] - altitude[i + 1])/scale_length);
    }
     return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
__global__ void overlap_parameter_d(int const num_layers, fp_t const * altitude,
                                    fp_t const scale_length, fp_t * alpha)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_layers - 1)
    {
        alpha[i] = exp(-1.*fabs(altitude[i] - altitude[i+1])/scale_length);
    }
    return;
}
#endif


/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
int sample_condensate(int const num_shape, int const num_x, fp_t const * x, fp_t const * y,
                      fp_t const * y_inverse, int const p, int const q, int const num_subcolumns,
                      int const num_layers, fp_t const * cloud_fraction,
                      fp_t const * lwc, fp_t const * iwc, fp_t const * overlap,
                      fp_t const * random1, fp_t const * random2, fp_t * ql, fp_t * qi)
{
    int i;
    for (i=0; i<num_subcolumns; ++i)
    {
        fp_t r[num_layers];
        cloudiness(num_layers, overlap, &(random1[i*num_layers]),
                   &(random2[i*(num_layers-1)]), r);
        int j;
        for (j=0; j<num_layers; ++j)
        {
            int offset = i*num_layers + j;
            if (r[j] > (1. - cloud_fraction[j]))
            {
                fp_t const qs = specific_saturation_humidity(num_shape, num_x,
                                                             x, y_inverse,
                                                             p, q, cloud_fraction[j]);
                fp_t const w = width(num_shape, num_x, x, y, p,
                                     q, cloud_fraction[j], lwc[j], iwc[j], qs);
                fp_t const b_inv = beta_inverse(num_shape, num_x, x,
                                                y_inverse, p, q, r[j]);
                fp_t const total_condensate = w*(b_inv - qs);
                fp_t const liquid_fraction = lwc[j]/(lwc[j] + iwc[j]);
                ql[offset] = total_condensate*liquid_fraction;
                qi[offset] = total_condensate*(1. - liquid_fraction);
            }
            else
            {
                ql[offset] = 0.;
                qi[offset] = 0.;
            }
        }
    }
    return GRTCODE_SUCCESS;
}


#ifdef __NVCC__
/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
__global__ void sample_condensate_d(int const num_shape, int const num_x, fp_t const * x,
                                    fp_t const * y, fp_t const * y_inverse, int const p,
                                    int const q, int const num_subcolumns,
                                    int const num_layers, fp_t const * cloud_fraction,
                                    fp_t const * lwc, fp_t const * iwc, fp_t const * overlap,
                                    fp_t const * random1, fp_t const * random2,
                                    fp_t * ql, fp_t * qi)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < num_subcolumns)
    {
        fp_t r[100];
        cloudiness(num_layers, overlap, &(random1[i*num_layers]),
                   &(random2[i*(num_layers-1)]), r);
        int j;
        for (j=0; j<num_layers; ++j)
        {
            int const offset = i*num_layers + j;
            if (r[j] > (1. - cloud_fraction[j]))
            {
                fp_t const qs = specific_saturation_humidity(num_shape, num_x,
                                                             x, y_inverse,
                                                             p, q, cloud_fraction[j]);
                fp_t const w = width(num_shape, num_x, x, y, p,
                                     q, cloud_fraction[j], lwc[j], iwc[j], qs);
                fp_t const b_inv = beta_inverse(num_shape, num_x, x,
                                                y_inverse, p, q, r[j]);
                fp_t const total_condensate = w*(b_inv - qs);
                fp_t const liquid_fraction = lwc[j]/(lwc[j] + iwc[j]);
                ql[offset] = total_condensate*liquid_fraction;
                qi[offset] = total_condensate*(1. - liquid_fraction);
            }
            else
            {
                ql[offset] = 0.;
                qi[offset] = 0.;
            }
        }
    }
    return;
}
#endif
