/* @brief Subcolumn generator for calculating stochastic clouds.*/
#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include "incomplete_beta.h"
#include "stochastic_clouds.h"


/* @brief Calculates random array for cloudiness from equation 1 from
          doi: 10.1256/qj.03.99.*/
static void cloudiness(int const num_layers, double const * overlap_parameter, double * x)
{
    int i;
    for (i=0; i<num_layers; ++i)
    {
        x[i] = ((double)(rand()))/((double)RAND_MAX);
    }
    double r[num_layers - 1];
    for (i=0; i<num_layers - 1; ++i)
    {
        r[i] = ((double)(rand()))/((double)RAND_MAX);
    }
    for (i=0; i<num_layers - 1; ++i)
    {
        if (r[i] <= overlap_parameter[i])
        {
            x[i + 1] = x[i];
        }
    }
    return;
}


/* @brief Calculates normalized saturation specific humidity from equation A1 from
          doi: 10.1175/MWR3257.1*/
static double specific_saturation_humidity(TotalWaterPDF_t const self,
                                           double const cloud_fraction)
{
    return beta_inverse(*(self.beta), self.p, self.q, 1. - cloud_fraction);
}


/* @brief Calculates the width of the total water probability distribution function (b - a)
          from equation A2 from doi: 10.1175/MWR3257.1, ignoring the parameter alpha.*/
static double width(TotalWaterPDF_t const self, double const cloud_fraction,
                    double const lwc, double const iwc, double const qs)
{
    return (lwc + iwc)/((((double)self.p)/((double)(self.p + self.q)))*
           (1. - beta_value(*(self.beta), self.p + 1, self.q, qs)) - qs*cloud_fraction);
}


/* @brief Constructs a TotalWaterPDF object.*/
void construct_water_pdf(TotalWaterPDF_t * self, int const p, int const q,
                         IncompleteBeta_t const * beta)
{
    self->p = p;
    self->q = q;
    self->beta = beta;
        fprintf(stderr, "p=%d q=%d self.p=%d self.q=%d\n",
            p, q, self->p, self->q);
    return;
}


/* @brief Destructs a TotalWaterPDF object.*/
void destruct_water_pdf(TotalWaterPDF_t * self)
{
    self->beta = NULL;
    return;
}


/* @brief Calculates the overlap parameter defined in equation 2 from
          doi: 10.1029/2004JD005100.*/
void overlap_parameter(int const num_layers, double const * altitude,
                       double const scale_length, double * alpha)
{
    int i;
    for (i=0; i<num_layers - 1; ++i)
    {
        alpha[i] = exp(-1.*fabs(altitude[i] - altitude[i + 1])/scale_length);
    }
    return;
}


/* @brief Draws samples of liquid and ice condensate mixing ratios from the total water
          (vapor + cloud liquid + ice) mixing ratio probability distribution function
          whose mean total condensate amount equals the sum of the input cloud liquid
          and ice condensate amounts and mean saturation humidity equals one minus
          the input cloud fraction for each layer.  This method is detailed in the
          appendix of doi: 10.1175/MWR3257.1.*/
void sample_condensate(TotalWaterPDF_t const self, int const num_layers,
                       double const * cloud_fraction, double const * lwc,
                       double const * iwc, double const * overlap, double * ql, double * qi)
{
    double x[num_layers];
    cloudiness(num_layers, overlap, x);
    int i;
    for (i=0; i<num_layers; ++i)
    {
        if (x[i] > (1. - cloud_fraction[i]))
        {
            double const qs = specific_saturation_humidity(self, cloud_fraction[i]);
            double const w = width(self, cloud_fraction[i], lwc[i], iwc[i], qs);
            double const total_condensate = w*(beta_inverse(*(self.beta), self.p, self.q,
                                                            x[i]) - qs);
            double const liquid_fraction = lwc[i]/(lwc[i] + iwc[i]);

            ql[i] = total_condensate*liquid_fraction;
            qi[i] = total_condensate*(1. - liquid_fraction);
        }
        else
        {
            ql[i] = 0.;
            qi[i] = 0.;
        }
    }
    return;
}