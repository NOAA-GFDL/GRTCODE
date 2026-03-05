#include <stdlib.h>
#include <string.h>
#include "optics_utils.h"
#include <stddef.h>  /* for NULL */
#include "cloud_pade_optics.h"

/* @brief Constructs an OpticalProperties object.*/
/* Constructs an OpticalProperties object by allocating two separate arrays for band_limits */
void construct_optics(OpticalProperties_t *self, int const num_bands, const double * const band_limits[2])

{
    fprintf(stdout, "construct_optics cloud optics with band [%d]\n", num_bands);
    self->num_bands = num_bands;

    /* Allocate a two-row array for band_limits */
     self->band_limits = malloc(2 * sizeof(double *));
    if (!self->band_limits) exit(1);  // Handle allocation failure

    self->band_limits[0] = malloc(num_bands * sizeof(double));  /* lower limits */
    self->band_limits[1] = malloc(num_bands * sizeof(double));  /* upper limits */
    if (!self->band_limits[0] || !self->band_limits[1]) exit(1);

    /* Copy from the provided 2D array */
    memcpy(self->band_limits[0], band_limits[0], num_bands * sizeof(double));
    memcpy(self->band_limits[1], band_limits[1], num_bands * sizeof(double));

    self->extinction_coefficient = malloc(num_bands * sizeof(double));
    self->single_scatter_albedo = malloc(num_bands * sizeof(double));
    self->asymmetry_factor = malloc(num_bands * sizeof(double));

    /* Optionally initialize bands array if needed,
       but currently it is not part of the structure per our header. */
}

/* Destructs an OpticalProperties object.*/
void destruct_optics(OpticalProperties_t *self)
{
    if(self->band_limits) {
        free(self->band_limits[0]);
        free(self->band_limits[1]);
        free(self->band_limits);
    }
    free(self->extinction_coefficient);
    free(self->single_scatter_albedo);
    free(self->asymmetry_factor);
}

/* 
 * Assume OpticalProperties_t is defined as follows:
 *
 * typedef struct {
 *     int num_bands;              // Number of coarse bands.
 *     fp_t *band_limits;          // Array of band limits of length (num_bands + 1).
 *     fp_t *extinction_coefficient;      // Coarse extinction coefficients (length: num_bands).
 *     fp_t *single_scatter_albedo;       // Coarse single-scatter albedo (length: num_bands).
 *     fp_t *asymmetry_factor;            // Coarse asymmetry factor (length: num_bands).
 * } OpticalProperties_t;
 */

/* Binary search helper: lower_bound
 * Returns the first index in arr[low..high] where arr[idx] >= target.
 */
static int lower_bound(const fp_t *arr, int low, int high, fp_t target) {
    int l = low, r = high + 1;
    while (l < r) {
        int mid = (l + r) / 2;
        if (arr[mid] < target)
            l = mid + 1;
        else
            r = mid;
    }
    return l;
}

/* Binary search helper: upper_bound
 * Returns the first index in arr[low..high] where arr[idx] > target.
 */
static int upper_bound(const fp_t *arr, int low, int high, fp_t target) {
    int l = low, r = high + 1;
    while (l < r) {
        int mid = (l + r) / 2;
        if (arr[mid] <= target)
            l = mid + 1;
        else
            r = mid;
    }
    return l-1;
}

/**
 * @brief Efficiently map coarse-band optical properties to a dense wavelength grid.
 *
 * This function assigns optical properties to output arrays over the entire
 * wavelength array. For each element:
 *   - If the wavelength is below the overall lower boundary of the coarse source,
 *     the first band’s optical properties are used.
 *   - If the wavelength is above the overall upper boundary, the last band’s properties are used.
 *   - For wavelengths between the coarse band's boundaries given by
 *     self.band_limits[input_ibnd] and self.band_limits[input_ibnd+1],
 *     the corresponding coarse-band properties are assigned.
 *
 * @param self       The coarse optical properties.
 * @param input_ibnd Pointer to the (0-indexed) coarse band index to use.
 *                   If NULL, defaults to 0.
 * @param wavenum    The dense (sorted) wavelength array (length N).
 * @param offset     The offset in the target arrays at which to begin assignments.
 * @param N          The number of wavelengths (length of wavenum and number of elements to assign).
 * @param beta       Output array to receive the extinction coefficients.
 * @param omega      Output array to receive the single-scatter albedo.
 * @param g          Output array to receive the asymmetry factor.
 */
void map_band_wave(OpticalProperties_t const self,
                   int const *input_ibnd,
                   const fp_t *wavenum,
                   int offset,
                   int N,
                   double *beta,
                   double *omega,
                   double *g)
{
    /* Determine the coarse-band index, defaulting to 0 if input_ibnd is NULL */
    int src_band = (input_ibnd == NULL) ? 0 : *input_ibnd;

    /* Overall boundaries from the coarse source:
       - Overall lower bound is self.band_limits[0]
       - Overall upper bound is self.band_limits[self.num_bands]
    */
    fp_t overall_low  = self.band_limits[0][0];
    fp_t overall_high = self.band_limits[1][self.num_bands-1];

    /* Coarse band boundaries for the selected band:
       - Lower bound is self.band_limits[src_band]
       - Upper bound is self.band_limits[src_band + 1]
    */
    fp_t band_low  = self.band_limits[0][src_band];
    fp_t band_high = self.band_limits[1][src_band];

    /* Use binary search to find the indices in the target dense grid.
       We search over the entire array (indices 0 through N-1) using the wavelength array.
    */
    int loIdx = lower_bound(wavenum, 0, N - 1, band_low);
    int hiIdx = upper_bound(wavenum, 0, N - 1, band_high);

    int j;

    /* Handle region before the first band */
    if (src_band == 0) {
        for (j = offset; j < offset + loIdx; j++) {
            beta[j]  = self.extinction_coefficient[0];
            omega[j] = self.single_scatter_albedo[0];
            g[j]     = self.asymmetry_factor[0];
        }
    }

    /* Main band-limited region */
    for (j = offset + loIdx; j < offset + hiIdx; j++) {
        beta[j]  = self.extinction_coefficient[src_band];
        omega[j] = self.single_scatter_albedo[src_band];
        g[j]     = self.asymmetry_factor[src_band];
    }

    /* Handle region beyond the last band */
    if (src_band == self.num_bands - 1) {
        for (j = offset + hiIdx; j < offset + N; j++) {
            beta[j]  = self.extinction_coefficient[src_band];
            omega[j] = self.single_scatter_albedo[src_band];
            g[j]     = self.asymmetry_factor[src_band];
        }
    }
}