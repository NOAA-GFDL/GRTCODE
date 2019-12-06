#ifndef DISORT_SHORTWAVE_H_
#define DISORT_SHORTWAVE_H_

#include "grtcode_utilities.h"


/** @brief Calculate upward and downward shortwave fluxes using DISORT.
    @return GRTCODE_SUCCESS or an error code.*/
EXTERN int disort_shortwave(Optics_t * const optics, /**< Optics.*/
                            fp_t const zen_dir, /**< Cosine of zenith angle.*/
                            fp_t * const surface_albedo, /**< Surface albedo (wavenumber).*/
                            fp_t const total_solar_irradiance, /**< Total solar irradiance [W m-2]. */
                            fp_t * const solar_flux, /**< Solar flux [cm] (wavenumber).*/
                            fp_t * const flux_up, /**< Upward flux [W cm m-2] (level, wavenumber).*/
                            fp_t * const flux_down /**< Downward flux [W cm m-2] (level, wavenumber).*/
                           );


#endif
