#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "netcdf.h"
#include "netcdf_utils.h"
#include <stddef.h>
#include "optics_utils.h"

#ifndef CLOUDS_LIB_H
#define CLOUDS_LIB_H

/**
 * @brief Structure for cloud optical properties.
 */
typedef struct {
    int m;            /**< Order of the polynomial */
    int n;            /**< Order of the polynomial */
    int nbnd;         /**< Number of bands */
    int nsizereg;     /**< Number of radius bins */

    fp_t rad_lwr;       /**< Radius lower bound [microns] */
    fp_t rad_upr;       /**< Radius upper bound [microns] */

    fp_t **band_lims_wvn; /**< Band limits (lower and upper bounds) */

    /* Pade approximant coefficients */
    fp_t ***pade_ext_p;
    fp_t ***pade_ssa_p;
    fp_t ***pade_asy_p;
    fp_t ***pade_ext_q;
    fp_t ***pade_ssa_q;
    fp_t ***pade_asy_q;

    /* Pade particle size regimes */
    fp_t **pade_sizereg;  /**< Pade radius bin limits [microns] */
    fp_t *pade_sizeref;   /**< Pade reference radii [microns] */
} ty_cloud_optics;

/**
 * @brief Initializes cloud optics from a NetCDF file.
 *
 * Reads spectral and size discretization, Pade coefficients, and band limits.
 *
 * @param cloud_optics Pointer to the cloud optics structure.
 * @param path Path to the NetCDF data file.
 * @return Status code (0 for success).
 */
int construct_cloud_optics(ty_cloud_optics *cloud_optics, const char *path);

/**
 * @brief Finalizes cloud optics structure.
 *
 * Frees allocated memory for cloud optics data.
 *
 * @param cloud_optics Pointer to the cloud optics structure.
 */
void finalize_ty_cloud_optics(ty_cloud_optics *cloud_optics);

/**
 * @brief Computes cloud optical properties using Pade approximants.
 *
 * @param cloud_optics Pointer to the cloud optics structure.
 * @param cwc Cloud water content.
 * @param re Effective radius.
 * @param optical_properties Output optical properties structure.
 * @param ibnd Band index to compute.
 * @return Status code (0 for success).
 */
int compute_all_from_pade(ty_cloud_optics *cloud_optics, fp_t cwc, fp_t re,
                          OpticalProperties_t *optical_properties, int ibnd);

/**
 * @brief Evaluates a Pade approximant.
 *
 * Computes numerator and denominator for a given set of coefficients.
 *
 * @param iband Band index.
 * @param nsizereg Number of radius bins.
 * @param n Order of polynomial numerator.
 * @param m Order of polynomial denominator.
 * @param irad Radius index.
 * @param re Effective radius.
 * @param pade_coeffs_p Numerator coefficients.
 * @param pade_coeffs_q Denominator coefficients.
 * @return Computed value of the Pade approximant.
 */
fp_t pade_eval_1(int iband, int nsizereg, int n, int m, int irad, fp_t re,
                 fp_t ***pade_coeffs_p, fp_t ***pade_coeffs_q);

#endif /* CLOUDS_LIB_H */
