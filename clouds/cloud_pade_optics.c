#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "cloud_pade_optics.h"
#include "netcdf.h"
#include "netcdf_utils.h"
#include "optics_utils.h"

// Define the ty_cloud_optics structure in the header file (do not redeclare it here)
#include "clouds_lib.h" // Ensure this header includes the correct ty_cloud_optics definition

// Function Prototypes
int construct_cloud_optics(ty_cloud_optics *self, const char *path);
void finalize_ty_cloud_optics(ty_cloud_optics *cloud_optics);
int compute_all_from_pade(ty_cloud_optics *cloud_optics, fp_t cwc, fp_t re, 
                          OpticalProperties_t *optical_properties, int ibnd);
fp_t pade_eval_1(int iband, int nsizereg, int n, int m, int irad, fp_t re, 
                 fp_t ***pade_coeffs_p, fp_t ***pade_coeffs_q);

// Function to initialize ty_cloud_optics
int construct_cloud_optics(ty_cloud_optics *self, const char *path) {
    if (!self || !path) {
        fprintf(stderr, "Error: null pointer passed to construct_cloud_optics.\n");
        return -1;
    }

    int ncid, retval;
    retval = nc_open(path, NC_NOWRITE, &ncid);
    if (retval != NC_NOERR) {
        fprintf(stderr, "NetCDF open failed: %s\n", nc_strerror(retval));
        return retval;
    }

    size_t dim_size;
    int dimid;

    // Read dimension sizes
    nc_inq_dimid(ncid, "Band", &dimid);
    nc_inq_dimlen(ncid, dimid, &dim_size);
    self->nbnd = (int)dim_size;

    nc_inq_dimid(ncid, "Re_range", &dimid);
    nc_inq_dimlen(ncid, dimid, &dim_size);
    self->nsizereg = (int)dim_size;

    nc_inq_dimid(ncid, "n", &dimid);
    nc_inq_dimlen(ncid, dimid, &dim_size);
    self->n = (int)dim_size;

    nc_inq_dimid(ncid, "m", &dimid);
    nc_inq_dimlen(ncid, dimid, &dim_size);
    self->m = (int)dim_size;

    // Allocate and read band limits
    self->band_lims_wvn = malloc(2 * sizeof(fp_t *));
    self->band_lims_wvn[0] = malloc(self->nbnd * sizeof(fp_t));
    self->band_lims_wvn[1] = malloc(self->nbnd * sizeof(fp_t));

    int varid;
    float *float_buffer = malloc(self->nbnd * sizeof(float));  // temp buffer

    nc_inq_varid(ncid, "Band_limits_lwr", &varid);
    nc_get_var_float(ncid, varid, float_buffer);
    for (int i = 0; i < self->nbnd; i++) self->band_lims_wvn[0][i] = (fp_t)float_buffer[i];

    nc_inq_varid(ncid, "Band_limits_upr", &varid);
    nc_get_var_float(ncid, varid, float_buffer);
    for (int i = 0; i < self->nbnd; i++) self->band_lims_wvn[1][i] = (fp_t)float_buffer[i];

    // Allocate and read radius regime
    self->pade_sizereg = malloc(2 * sizeof(fp_t *));
    self->pade_sizereg[0] = malloc(self->nsizereg * sizeof(fp_t));
    self->pade_sizereg[1] = malloc(self->nsizereg * sizeof(fp_t));
    self->pade_sizeref   = malloc(self->nsizereg * sizeof(fp_t));

    nc_inq_varid(ncid, "Effective_Radius_limits_lwr", &varid);
    nc_get_var_float(ncid, varid, float_buffer);
    for (int i = 0; i < self->nsizereg; i++) self->pade_sizereg[0][i] = (fp_t)float_buffer[i];

    nc_inq_varid(ncid, "Effective_Radius_limits_upr", &varid);
    nc_get_var_float(ncid, varid, float_buffer);
    for (int i = 0; i < self->nsizereg; i++) self->pade_sizereg[1][i] = (fp_t)float_buffer[i];

    nc_inq_varid(ncid, "Effective_Radius_Ref", &varid);
    nc_get_var_float(ncid, varid, float_buffer);
    for (int i = 0; i < self->nsizereg; i++) self->pade_sizeref[i] = (fp_t)float_buffer[i];

    // Allocate 3D arrays for Pade coefficients
    size_t num = self->nbnd;
    size_t nr = self->nsizereg;
    size_t n  = self->n;
    size_t m  = self->m;

    #define ALLOC_P3(name) \
        self->name = malloc(num * sizeof(fp_t **)); \
        for (int ib = 0; ib < num; ++ib) { \
            self->name[ib] = malloc(nr * sizeof(fp_t *)); \
            for (int ir = 0; ir < nr; ++ir) \
                self->name[ib][ir] = malloc((strcmp(#name, "pade_ext_q") == 0 || strcmp(#name, "pade_ssa_q") == 0 || strcmp(#name, "pade_asy_q") == 0 ? m : n) * sizeof(fp_t)); \
        }

    ALLOC_P3(pade_ext_p)
    ALLOC_P3(pade_ext_q)
    ALLOC_P3(pade_ssa_p)
    ALLOC_P3(pade_ssa_q)
    ALLOC_P3(pade_asy_p)
    ALLOC_P3(pade_asy_q)

    #define READ_P3(nc_varname, target, order) \
    nc_inq_varid(ncid, nc_varname, &varid); \
    nc_get_var_float(ncid, varid, p3_data); \
    for (int ib = 0; ib < num; ++ib) \
        for (int ir = 0; ir < nr; ++ir) \
            for (int icoef = 0; icoef < order; ++icoef) { \
                int idx = ib + ir * num + icoef * nr * num; \
                target[ib][ir][icoef] = (fp_t)p3_data[idx]; \
            }

    float *p3_data = malloc(self->n * self->nsizereg * self->nbnd * sizeof(float));  // safe upper bound
    READ_P3("Pade_ext_p", self->pade_ext_p, n);
    READ_P3("Pade_ext_q", self->pade_ext_q, m);
    READ_P3("Pade_ssa_p", self->pade_ssa_p, n);
    READ_P3("Pade_ssa_q", self->pade_ssa_q, m);
    READ_P3("Pade_asy_p", self->pade_asy_p, n);
    READ_P3("Pade_asy_q", self->pade_asy_q, m);

    // Set outer radius bounds
    self->rad_lwr = self->pade_sizereg[0][0];
    self->rad_upr = self->pade_sizereg[1][self->nsizereg - 1];

    free(float_buffer);
    nc_close(ncid);
    return 0;
}


// Function to finalize ty_cloud_optics
void finalize_ty_cloud_optics(ty_cloud_optics *cloud_optics) {
    if (!cloud_optics) return;

    free(cloud_optics->pade_ext_p);
    free(cloud_optics->pade_ssa_p);
    free(cloud_optics->pade_asy_p);
    free(cloud_optics->pade_ext_q);
    free(cloud_optics->pade_ssa_q);
    free(cloud_optics->pade_asy_q);

    free(cloud_optics->pade_sizereg[0]);
    free(cloud_optics->pade_sizereg[1]);
    free(cloud_optics->pade_sizeref);
    free(cloud_optics->band_lims_wvn[0]);
    free(cloud_optics->band_lims_wvn[1]);
    free(cloud_optics->band_lims_wvn);
}

// Function to compute optical properties from Pade coefficients
int compute_all_from_pade(ty_cloud_optics *cloud_optics, fp_t cwc, fp_t re, 
                          OpticalProperties_t *optical_properties, int ibnd) {


    if (cwc > 0) {
        fp_t radius = re;
        int irad = 0;
        for (; irad < cloud_optics->nsizereg; ++irad) {
            if (cloud_optics->pade_sizereg[0][irad] <= radius &&
                cloud_optics->pade_sizereg[1][irad] >= radius) {
                break;
            }
        }
        if (irad >= cloud_optics->nsizereg) {
            optical_properties->extinction_coefficient[ibnd] = 0.0;
            optical_properties->single_scatter_albedo[ibnd] = 0.0;
            optical_properties->asymmetry_factor[ibnd] = 0.0;
            return 0;
        }

        fp_t re_offset = radius - cloud_optics->pade_sizeref[irad];

        optical_properties->extinction_coefficient[ibnd] = cwc * pade_eval_1(
            ibnd, cloud_optics->nsizereg, cloud_optics->n, cloud_optics->m, irad,
            re_offset, cloud_optics->pade_ext_p, cloud_optics->pade_ext_q);

        optical_properties->single_scatter_albedo[ibnd] = pade_eval_1(
            ibnd, cloud_optics->nsizereg, cloud_optics->n, cloud_optics->m, irad,
            re_offset, cloud_optics->pade_ssa_p, cloud_optics->pade_ssa_q);

        optical_properties->asymmetry_factor[ibnd] = pade_eval_1(
            ibnd, cloud_optics->nsizereg, cloud_optics->n, cloud_optics->m, irad,
            re_offset, cloud_optics->pade_asy_p, cloud_optics->pade_asy_q);
    } else {
        optical_properties->extinction_coefficient[ibnd] = 0.0;
        optical_properties->single_scatter_albedo[ibnd] = 0.0;
        optical_properties->asymmetry_factor[ibnd] = 0.0;
    }


    return 0;
}

// Implementation in pade_eval.c
fp_t pade_eval_1(int iband, int nsizereg, int n, int m, int irad, fp_t re, 
                 fp_t ***pade_coeffs_p, fp_t ***pade_coeffs_q) {
    fp_t numer = pade_coeffs_p[iband][irad][0];
    for (int i = 1; i < n; ++i) {
        numer = pade_coeffs_p[iband][irad][i] + re * numer;
    }
    fp_t denom = pade_coeffs_q[iband][irad][0];
    for (int i = 1; i < m; ++i) {
        denom = pade_coeffs_q[iband][irad][i] + re * denom;
    }
    return numer / denom;
}
