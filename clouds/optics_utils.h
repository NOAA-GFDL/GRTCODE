#ifndef OPTICS_UTILS_H
#define OPTICS_UTILS_H

#include <stddef.h>  /* for NULL */

/**
 * @brief OpticalProperties_t
 *
 * This structure represents the coarse-band optical properties.
 * In this design, the field `band_limits` is expected to be a two‐row array:
 *   - band_limits[0]: Lower limits of each coarse band (length: num_bands)
 *   - band_limits[1]: Upper limits of each coarse band (length: num_bands)
 *
 * The overall lower bound is band_limits[0][0] and the overall upper bound is
 * band_limits[1][num_bands-1].
 *
 * The optical properties arrays (extinction_coefficient, single_scatter_albedo,
 * and asymmetry_factor) contain one value per coarse band.
 */
typedef struct {
    int num_bands;                /* Number of coarse bands */
    double **band_limits;         /* 2D array of band limits (2 x num_bands) */
    double *extinction_coefficient;   /* Coarse extinction coefficients (length: num_bands) */
    double *single_scatter_albedo;    /* Coarse single-scatter albedo (length: num_bands) */
    double *asymmetry_factor;         /* Coarse asymmetry factor (length: num_bands) */
    /* Optionally, additional members can be added here */
} OpticalProperties_t;

/**
 * @brief Constructs an OpticalProperties object.
 *
 * Initializes an OpticalProperties_t instance with a given number of coarse bands.
 * The source array 'band_limits_src' should contain the band limits arranged so that
 * the first num_bands entries are the lower limits and the next num_bands entries are
 * the corresponding upper limits.
 *
 * @param self             Pointer to the OpticalProperties_t instance.
 * @param num_bands        Number of coarse bands.
 * @param band_limits      Source array containing the band limits (length: 2 * num_bands).
 */
void construct_optics(OpticalProperties_t *self, int num_bands, const double *const band_limits[2]);

/**
 * @brief Destructs an OpticalProperties object.
 *
 * Frees any memory allocated within the OpticalProperties_t.
 *
 * @param self  Pointer to the OpticalProperties_t instance.
 */
void destruct_optics(OpticalProperties_t *self);

/**
 * @brief Efficiently map coarse-band optical properties onto a dense wavelength grid.
 *
 * Given a sorted wavelength array (wavenum) of length N and output arrays (beta,
 * omega, and g) that start at a given offset, this function fills in the optical
 * properties as follows:
 *
 *   - For wavelengths below the overall lower bound (band_limits[0][0]), the
 *     properties from band 0 are assigned.
 *
 *   - For wavelengths above the overall upper bound (band_limits[1][num_bands-1]),
 *     the properties from the last band are assigned.
 *
 *   - For wavelengths within a given coarse band's limits—that is, between
 *     band_limits[0][input_ibnd] and band_limits[1][input_ibnd]—the corresponding
 *     coarse-band properties are assigned.
 *
 * @param self       The coarse optical properties.
 * @param input_ibnd Pointer to the (0-indexed) coarse band index to use.
 *                   If NULL, defaults to 0.
 * @param wavenum    The sorted dense wavelength array (length N).
 * @param offset     The offset in the output arrays at which assignments begin.
 * @param N          The number of wavelengths (length of wavenum and number of values to assign).
 * @param beta       Output array for extinction coefficients.
 * @param omega      Output array for single-scatter albedo.
 * @param g          Output array for asymmetry factor.
 */
void map_band_wave(OpticalProperties_t const self,
                   int const *input_ibnd,
                   const double *wavenum,
                   int offset,
                   int N,
                   double *beta,
                   double *omega,
                   double *g);

#endif /* OPTICS_UTILS_H */
