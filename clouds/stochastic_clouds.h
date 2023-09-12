#ifndef STOCHASTIC_CLOUDS_H
#define STOCHASTIC_CLOUDS_H


#include "debug.h"
#include "incomplete_beta.h"


/* @brief Total water probability density function.*/
typedef struct TotalWaterPDF
{
    Device_t device; /*Device to run the code on.*/
    int num_shape; /*Number of possible beta distribution shape parameters.*/
    int num_x; /*Number of beta distribution table values per shape parameter.*/
    fp_t * x; /*Table input values.*/
    fp_t * y; /*Table of calculated values.*/
    fp_t * y_inverse; /*Table of inverse values.*/
    int p; /*Incomplete beta distribution shape parameter.*/
    int q; /*Incomplete beta distribution shape parameter.*/
} TotalWaterPDF_t;


/* @brief Constructs a TotalWaterPDF object.*/
EXTERN int create_water_pdf(
    TotalWaterPDF_t * self, /**< Total water PDF object.*/
    int const p, /**< Beta distribution shape parameter.*/
    int const q, /**< Beta distribution shape parameter.*/
    IncompleteBeta_t const * beta /**< Beta distribution object.*/
);


/* @brief Destructs a TotalWaterPDF object.*/
EXTERN int destroy_water_pdf(
    TotalWaterPDF_t * self /**< Total water PDF object.*/
);


/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
int overlap_parameter(
    int const num_layers, /**< Number of layers.*/
    fp_t const * altitude, /**< Altitude [km] (layer).*/
    fp_t const scale_length, /**< Scale length [km].*/
    fp_t * alpha /**< Overlap parameter (layer - 1).*/
);


#ifdef __NVCC__
__global__ void overlap_parameter_d(
    int const num_layers, /**< Number of layers.*/
    fp_t const * altitude, /**< Altitude [km] (layer).*/
    fp_t const scale_length, /**< Scale length [km].*/
    fp_t * alpha /**< Overlap parameter (layer - 1).*/
);
#endif


/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
int sample_condensate(
    int const num_shape, /**< Number of beta distribution shape parameter.*/
    int const num_x, /**< Number of beta distribution x values.*/
    fp_t const * x, /**< Beta distribution x values (x).*/
    fp_t const * y, /**< Beta distribution values (shape, shape, x).*/
    fp_t const * y_inverse, /**< Beta distribution deviates (shape, shape, x).*/
    int const p, /**< Beta distribution shape parameter.*/
    int const q, /**< Beta distribution shape parameter.*/
    int const num_subcolumns, /**< Number of subcolumns.*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * cloud_fraction, /**< Grid-cell mean cloud fraction (layer).*/
    fp_t const * lwc, /**< Grid-cell mean liquid condensate (layer).*/
    fp_t const * iwc, /**< Grid-cell mean ice condensate (layer).*/
    fp_t const * overlap, /**< Overlap parameter (layer - 1).*/
    fp_t const * random1, /**< Random numbers (subcolumn, layer).*/
    fp_t const * random2, /**< Random numbers (subcolumn, layer - 1).*/
    fp_t * ql, /**< Liquid condensate (subcolumn, layer).*/
    fp_t * qi /**< Ice condensate (subcolumn, layer).*/
);


#ifdef __NVCC__
/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
__global__ void sample_condensate_d(
    int const num_shape, /**< Number of beta distribution shape parameter.*/
    int const num_x, /**< Number of beta distribution x values.*/
    fp_t const * x, /**< Beta distribution x values (x).*/
    fp_t const * y, /**< Beta distribution values (shape, shape, x).*/
    fp_t const * y_inverse, /**< Beta distribution deviates (shape, shape, x).*/
    int const p, /**< Beta distribution shape parameter.*/
    int const q, /**< Beta distribution shape parameter.*/
    int const num_subcolumns, /**< Number of subcolumns.*/
    int const num_layers, /**< Number of layers.*/
    fp_t const * cloud_fraction, /**< Grid-cell mean cloud fraction (layer).*/
    fp_t const * lwc, /**< Grid-cell mean liquid condensate (layer).*/
    fp_t const * iwc, /**< Grid-cell mean ice condensate (layer).*/
    fp_t const * overlap, /**< Overlap parameter (layer - 1).*/
    fp_t const * random1, /**< Random numbers (subcolumn, layer).*/
    fp_t const * random2, /**< Random numbers (subcolumn, layer - 1).*/
    fp_t * ql, /**< Liquid condensate (subcolumn, layer).*/
    fp_t * qi /**< Ice condensate (subcolumn, layer).*/
);
#endif


#endif
