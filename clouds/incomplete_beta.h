#ifndef INCOMPLETE_BETA_H
#define INCOMPLETE_BETA_H

#include "debug.h"
#include "device.h"


/* @brief Incomplete beta distribution look-up tables.*/
typedef struct IncompleteBeta
{
    Device_t device; /*Device to run the code on.*/
    int num_shape; /*Number of possible beta distribution shape parameters.*/
    int num_x; /*Number of beta distribution table values per shape parameter.*/
    int * p; /*Incomplete beta distribution shape parameters.*/
    int * q; /*Incomplete beta distribution shape parameters.*/
    fp_t * x; /*Table input values.*/
    fp_t * y; /*Table of calculated values.*/
    fp_t * y_inverse; /*Table of inverse values.*/
} IncompleteBeta_t;


/* @brief Constructs an IncompleteBeta object.*/
EXTERN int create_beta(IncompleteBeta_t * self,
                       char const * path,
                       Device_t const * const device);


/* @brief Destructs an IncompleteBeta object.*/
EXTERN int destroy_beta(IncompleteBeta_t * self);


/* @brief Calculates the inverse of the incomplete beta function.*/
HOST DEVICE fp_t beta_inverse(int const num_shape,
                              int const num_x,
                              fp_t const * x,
                              fp_t const * y_inverse,
                              int const p,
                              int const q,
                              fp_t const z);


/* @brief Calculates incomplete beta function.*/
HOST DEVICE fp_t beta_value(int const num_shape,
                            int const num_x,
                            fp_t const * x,
                            fp_t const * y,
                            int const p,
                            int const q,
                            fp_t const z);


#endif
