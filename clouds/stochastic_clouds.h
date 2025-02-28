#ifndef STOCHASTIC_CLOUDS_H
#define STOCHASTIC_CLOUDS_H

#include "incomplete_beta.h"


/* @brief Total water probability density function.*/
typedef struct TotalWaterPDF
{
    IncompleteBeta_t const * beta; /*Incomplete beta distribution.*/
    int p; /*Incomplete beta distribution shape parameter.*/
    int q; /*Incomplete beta distribution shape parameter.*/
} TotalWaterPDF_t;


/* @brief Constructs a TotalWaterPDF object.*/
void construct_water_pdf(TotalWaterPDF_t * self, int const p, int const q,
                         IncompleteBeta_t const * beta);


/* @brief Destructs a TotalWaterPDF object.*/
void destruct_water_pdf(TotalWaterPDF_t * self);


/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
void overlap_parameter(int const num_layers, double const * altitude,
                       double const scale_length, double * alpha);


/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
void sample_condensate(TotalWaterPDF_t const self, int const num_layers,
                       double const * cloud_fraction, double const * lwc,
                       double const * iwc, double const * overlap, double * ql, double * qi);


#endif
