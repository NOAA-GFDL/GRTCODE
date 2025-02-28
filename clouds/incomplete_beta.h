#ifndef INCOMPLETE_BETA_H
#define INCOMPLETE_BETA_H


/* @brief Incomplete beta distribution look-up tables.*/
typedef struct IncompleteBeta
{
    int num_shape; /*Number of possible beta distribution shape parameters.*/
    int num_x; /*Number of beta distribution table values per shape parameter.*/
    int * p; /*Incomplete beta distribution shape parameters.*/
    int * q; /*Incomplete beta distribution shape parameters.*/
    double * x; /*Table input values.*/
    double * y; /*Table of calculated values.*/
    double * y_inverse; /*Table of inverse values.*/
} IncompleteBeta_t;


/* @brief Constructs an IncompleteBeta object.*/
void construct_beta(IncompleteBeta_t * self, char const * path);


/* @brief Destructs an IncompleteBeta object.*/
void destruct_beta(IncompleteBeta_t * self);


/* @brief Calculates the inverse of the incomplete beta function.*/
double beta_inverse(IncompleteBeta_t const self, int const p, int const q, double const x);


/* @brief Calculates incomplete beta function.*/
double beta_value(IncompleteBeta_t const self, int const p, int const q, double const x);


#endif
