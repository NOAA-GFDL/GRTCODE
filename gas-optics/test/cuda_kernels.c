/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifdef __NVCC__

#include <math.h>
#include <stdint.h>
#include "cuda_kernels.cuh"
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "kernel_utils.h"
#include "line_shape.h"
#include "RFM_voigt.h"
#include "spectral_bin.h"
#include "spectral_bin-internal.h"
#include "tips2017.h"


#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ < 600
/*Atomic add for doubles must be defined for some older GPUs.*/
__device__ static double atomicAdd(double *address, /*Address to add value to.*/
                                   double val /*Value to add.*/
                                  )
{
    unsigned long long int *address_as_ull = (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull;
    unsigned long long int assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
#endif


/*Calculate pressure-shifted line center positions.*/
__global__ void calc_line_centers_d(uint64_t const num_lines, int const num_layers,
                                    fp_t const * const v0, fp_t const * const delta,
                                    fp_t const * const p, fp_t * const vnn)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        int i;
        for (i=0;i<num_layers;++i)
        {
            vnn[i*num_lines+j] = v0[j] + delta[j]*p[i];
        }
    }
    return;
}


/*Calculate the total partition functions for each isotopologue.*/
__global__ void calc_partition_functions_d(int const num_layers, int const mol_id,
                                           int const num_iso, fp_t const * const t,
                                           fp_t * const q)
{
    int const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_iso)
    {
        int i;
        for (i=0;i<num_layers;++i)
        {
            q[i*num_iso+j] = 1.f/Q(mol_id,t[i],j+1);
        }
    }
    return;
}


/*Calculate temperature-corrected line intensities.*/
__global__ void calc_line_strengths_d(uint64_t const num_lines, int const num_layers,
                                      int const num_iso, int const * const iso,
                                      fp_t const * const s0, fp_t const * const vnn,
                                      fp_t const * const en, fp_t const * const t,
                                      fp_t const * const q,  fp_t * const snn)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        fp_t const c2 = -1.4387686f;
        int i;
        for (i=0;i<num_layers;++i)
        {
            snn[i*num_lines+j] = s0[j]*EXP(c2*en[j]/t[i])*
                                 (1.f - EXP(c2*vnn[j]/t[i]))*
                                 q[i*num_iso+iso[j]-1];
        }
    }
    return;
}


/*Calculate lorentz halfwidths.*/
__global__ void calc_lorentz_hw_d(uint64_t const num_lines, int const num_layers,
                                  fp_t const * const n, fp_t const * const yair,
                                  fp_t const * const yself, fp_t const * const t,
                                  fp_t const * const p, fp_t const * const ps,
                                  fp_t * const gamma)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        fp_t const tref = 296.f;
        int i;
        for (i=0;i<num_layers;++i)
        {
            gamma[i*num_lines+j] = POW(tref/t[i],n[j])*(yair[j]*(p[i]-ps[i]) +
                                   yself[j]*ps[i]);
        }
    }
    return;
}


/*Calculate doppler halfwidths.*/
__global__ void calc_doppler_hw_d(uint64_t const num_lines, int const num_layers,
                                  fp_t const m, fp_t const * const vnn,
                                  fp_t const * const t, fp_t * const alpha)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        fp_t const sqrt_ln2 = 0.83255461115f;
        fp_t const kb = 1.380658E-16; /*[erg/K]*/
        fp_t const c = 2.99792458E10; /*[cm/s]*/
        int i;
        for (i=0;i<num_layers;++i)
        {
            alpha[i*num_lines+j] = sqrt_ln2*vnn[i*num_lines+j]*
                                   SQRT((2.f*kb*t[i])/(m*c*c));
        }
    }
    return;
}


/*Sort the line parameters in order of line center wavenumber.*/
__global__ void sort_lines_d(uint64_t const num_lines, int const num_layers,
                             fp_t * const vnn, fp_t * const snn,
                             fp_t * const gamma, fp_t * const alpha)
{
    int const k = blockIdx.x*blockDim.x + threadIdx.x;
    if (k < num_layers)
    {
        fp_t *v = &(vnn[k*num_lines]);
        fp_t *s = &(snn[k*num_lines]);
        fp_t *g = &(gamma[k*num_lines]);
        fp_t *a = &(alpha[k*num_lines]);
        uint64_t i;
        for (i=1;i<num_lines;++i)
        {
            fp_t value[4];
            value[0] = v[i];
            value[1] = s[i];
            value[2] = g[i];
            value[3] = a[i];
            uint64_t j = i;
            while (j > 0 && v[j-1] > value[0])
            {
                v[j] = v[j-1];
                s[j] = s[j-1];
                g[j] = g[j-1];
                a[j] = a[j-1];
                j--;
            }
            if (i != j)
            {
                v[j] = value[0];
                s[j] = value[1];
                g[j] = value[2];
                a[j] = value[3];
            }
        }
    }
    return;
}


/*Calculate optical depths.*/
__global__ void calc_optical_depth_bin_sweep_d(uint64_t const num_lines, int const num_layers,
                                               fp_t * const vnn, fp_t * const snn,
                                               fp_t * const gamma, fp_t * const alpha,
                                               fp_t const * const n, SpectralBins_t bins,
                                               fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < bins.n)
    {
        int i;
        for (i=0;i<num_layers;++i)
        {
            fp_t const *v = &(vnn[i*num_lines]);
            fp_t const *s = &(snn[i*num_lines]);
            fp_t const *g = &(gamma[i*num_lines]);
            fp_t const *a = &(alpha[i*num_lines]);
            uint64_t nbin_local = 1;
            uint64_t nbin_remote = 25;

            /*Find the "local" lines.*/
            uint64_t nbin = nbin_local;
            fp_t leftw = j > nbin ? bins.w[NIP*(j-nbin)] : bins.w[0];
            fp_t rightw = j >= (bins.n-1)-nbin ? bins.w[NIP*bins.n-1]
                          : bins.w[NIP*(j+nbin+1)-1];

            uint64_t left;
            uint64_t right;
            if (leftw <= v[num_lines-1] && rightw >= v[0])
            {
                uint64_t tmp;
                bracket(num_lines,
                        v,
                        leftw,
                        &left,
                        &tmp);
                bracket(num_lines - left,
                        &(v[left]),
                        rightw,
                        &tmp,
                        &right);
                right += left;

                LineShapeInputs_t in;
                in.wres = bins.wres;
                uint64_t const blocksize = 16;
                fp_t t[blocksize];
                uint64_t const num_blocks = (bins.r[j]-bins.l[j]+1)/blocksize;
                uint64_t const remainder = (bins.r[j]-bins.l[j]+1)%blocksize;
                uint64_t k;
                for (k=left;k<=right;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    in.num_wpoints = blocksize;
                    uint64_t b;
                    for (b=0;b<num_blocks;++b)
                    {
                        in.w = bins.w0 + (bins.l[j] + b*blocksize)*bins.wres;
                        rfm_voigt_line_shape(in,
                                             t);
                        uint64_t l;
                        for (l=0;l<blocksize;l++)
                        {
                            atomicAdd(&(tau[i*bins.num_wpoints+bins.l[j]+b*blocksize+l]),
                                      s[k]*n[i]*t[l]);
                        }
                    }
                    if (remainder > 0)
                    {
                        in.w = bins.w0 + (bins.l[j] + num_blocks*blocksize)*bins.wres;
                        in.num_wpoints = remainder;
                        rfm_voigt_line_shape(in,
                                             t);
                        uint64_t l;
                        for (l=0;l<remainder;l++)
                        {
                            atomicAdd(&(tau[i*bins.num_wpoints+bins.l[j]+num_blocks*blocksize+l]),
                                      s[k]*n[i]*t[l]);
                        }
                    }
                }
            }
            else if (leftw > v[num_lines-1])
            {
                left = num_lines;
            }
            else
            {
                right = (uint64_t)(-1);
            }

            /*Find the "remote" lines.*/
            nbin = nbin_remote;
            fp_t leftw_r = j > nbin ? bins.w[NIP*(j-nbin)] : bins.w[0];
            if (leftw >= v[0] && leftw_r <= v[num_lines-1])
            {
                uint64_t left_r;
                uint64_t tmp;
                bracket(left,
                        v,
                        leftw_r,
                        &left_r,
                        &tmp);
                LineShapeInputs_t in;
                in.w = bins.w[j*NIP];
                in.num_wpoints = NIP;
                in.wres = bins.w[j*NIP+1] - in.w;
                fp_t t_r[NIP];
                uint64_t k;
                for (k=left_r;k<left;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins.n*NIP + j*NIP + l;
                        bins.tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }

            fp_t rightw_r = j >= (bins.n-1)-nbin ? bins.w[NIP*bins.n-1]
                            : bins.w[NIP*(j+nbin+1)-1];
            if (rightw <= v[num_lines-1] && rightw_r >= v[0])
            {
                uint64_t f = 0;
                if (right == (uint64_t)(-1))
                {
                    f = 1;
                }
                uint64_t right_r;
                uint64_t tmp;
                bracket(num_lines - (right + f),
                        &(v[right+f]),
                        rightw_r,
                        &tmp,
                        &right_r);
                right_r += right + f;
                LineShapeInputs_t in;
                in.w = bins.w[j*NIP];
                in.num_wpoints = NIP;
                in.wres = bins.w[j*NIP+1] - in.w;
                fp_t t_r[NIP];
                uint64_t k;
                for (k=right+1;k<=right_r;++k)
                {
                    in.line_center = v[k];
                    in.lorentz_hwhm = g[k];
                    in.doppler_hwhm = a[k];
                    rfm_voigt_line_shape(in,
                                         t_r);
                    uint64_t l;
                    for (l=0;l<NIP;++l)
                    {
                        uint64_t offset = i*bins.n*NIP + j*NIP + l;
                        bins.tau[offset] += s[k]*n[i]*t_r[l];
                    }
                }
            }
        }
    }
    return;
}


/*Calculate optical depths.*/
__global__ void calc_optical_depth_line_sweep_d(uint64_t const num_lines, int const num_layers,
                                                fp_t * const vnn, fp_t * const snn,
                                                fp_t * const gamma, fp_t * const alpha,
                                                fp_t const * const n, SpectralBins_t bins,
                                                fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        fp_t const bin_width = bins.wres*bins.ppb;
        int i;
        for (i=0;i<num_layers;++i)
        {
            uint64_t o = i*num_lines + j;
            LineShapeInputs_t in;
            in.line_center = vnn[o];
            in.lorentz_hwhm = gamma[o];
            in.doppler_hwhm = alpha[o];
            in.wres = bins.wres;

            /*Local lines.*/
            fp_t wcutoff = 1.5f;
            fp_t leftw = in.line_center - wcutoff;
            if (leftw < bins.w0)
            {
                leftw = bins.w0;
            }
            uint64_t left = floor((leftw-bins.w0)/bin_width);

            fp_t rightw = in.line_center + wcutoff;
            fp_t maxw = bins.w0 + bins.num_wpoints*bins.wres;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right = floor((rightw-bins.w0)/bin_width);

            uint64_t k;
            for (k=left;k<=right;++k)
            {
                uint64_t const blocksize = 16;
                fp_t t[blocksize];
                uint64_t const num_blocks = (bins.r[k]-bins.l[k]+1)/blocksize;
                uint64_t const remainder = (bins.r[k]-bins.l[k]+1)%blocksize;
                in.num_wpoints = blocksize;
                uint64_t b;
                for (b=0;b<num_blocks;++b)
                {
                    in.w = bins.w0 + (bins.l[k] + b*blocksize)*bins.wres;
                    rfm_voigt_line_shape(in,
                                         t);
                    uint64_t l;
                    for (l=0;l<blocksize;l++)
                    {
                        atomicAdd(&(tau[i*bins.num_wpoints+bins.l[k]+b*blocksize+l]),
                                  snn[o]*n[i]*t[l]);
                    }
                }
                if (remainder > 0)
                {
                    in.w = bins.w0 + (bins.l[k] + num_blocks*blocksize)*bins.wres;
                    in.num_wpoints = remainder;
                    rfm_voigt_line_shape(in,
                                         t);
                    uint64_t l;
                    for (l=0;l<remainder;l++)
                    {
                        atomicAdd(&(tau[i*bins.num_wpoints+bins.l[k]+num_blocks*blocksize+l]),
                                  snn[o]*n[i]*t[l]);
                    }
                }
            }

            /*Remote lines.*/
            wcutoff = 25.f;
            leftw = in.line_center - wcutoff;
            if (leftw < bins.w0)
            {
                leftw = bins.w0;
            }
            uint64_t left_r = floor((leftw-bins.w0)/bin_width);
            in.num_wpoints = NIP;
            fp_t t_r[NIP];
            for (k=left_r;k<left;++k)
            {
                in.w = bins.w[k*NIP];
                in.wres = bins.w[k*NIP+1] - in.w;
                rfm_voigt_line_shape(in,
                                     t_r);
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins.n*NIP + k*NIP + l;
                    atomicAdd(&(bins.tau[offset]),
                              snn[o]*n[i]*t_r[l]);
                }
            }

            rightw = in.line_center + wcutoff;
            if (rightw > maxw)
            {
                rightw = maxw;
            }
            uint64_t right_r = floor((rightw-bins.w0)/bin_width);
            for (k=right+1;k<=right_r;++k)
            {
                in.w = bins.w[k*NIP];
                in.wres = bins.w[k*NIP+1] - in.w;
                rfm_voigt_line_shape(in,
                                     t_r);
                uint64_t l;
                for (l=0;l<NIP;++l)
                {
                    uint64_t offset = i*bins.n*NIP + k*NIP + l;
                    atomicAdd(&(bins.tau[offset]),
                              snn[o]*n[i]*t_r[l]);
                }
            }
        }
    }
    return;
}


/*Calculate optical depths by sampling all lines.*/
__global__ void calc_optical_depth_line_sample_d(uint64_t const num_lines, int const num_layers,
                                                 fp_t * const vnn, fp_t * const snn,
                                                 fp_t * const gamma, fp_t * const alpha,
                                                 fp_t const * const n, SpectralBins_t const bins,
                                                 fp_t * const tau)
{
    uint64_t j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_lines)
    {
        uint64_t const fsteps = ceil(25.f/bins.wres);
        int i;
        for (i=0;i<num_layers;++i)
        {
            uint64_t loffset = i*num_lines + j;
            LineShapeInputs_t in;
            in.line_center = vnn[loffset];
            in.lorentz_hwhm = gamma[loffset];
            in.doppler_hwhm = alpha[loffset];
            in.wres = bins.wres;
            uint64_t fcenterid = floor((2*((in.line_center-bins.w0)/
                                       bins.wres)+1)/2);
            if (fcenterid < bins.num_wpoints)
            {
                uint64_t s = (int64_t)(fcenterid-fsteps) < 0 ? 0 :
                             fcenterid-fsteps;
                uint64_t e = fcenterid+fsteps >= bins.num_wpoints ?
                             bins.num_wpoints-1 : fcenterid+fsteps;
                uint64_t const blocksize = 16;
                fp_t t[blocksize];
                uint64_t const num_blocks = (e-s+1)/blocksize;
                uint64_t const remainder = (e-s+1)%blocksize;
                in.num_wpoints = blocksize;
                uint64_t k;
                for (k=0;k<num_blocks;++k)
                {
                    in.w = (s + k*blocksize)*bins.wres + bins.w0;
                    rfm_voigt_line_shape(in,
                                         t);
                    uint64_t l;
                    for (l=0;l<blocksize;l++)
                    {
                        atomicAdd(&(tau[i*bins.num_wpoints+s+k*blocksize+l]),
                                  snn[loffset]*n[i]*t[l]);
                    }
                }
                if (remainder > 0)
                {
                    in.w = (s + num_blocks*blocksize)*bins.wres + bins.w0;
                    in.num_wpoints = remainder;
                    rfm_voigt_line_shape(in,
                                         t);
                    uint64_t l;
                    for (l=0;l<remainder;l++)
                    {
                        atomicAdd(&(tau[i*bins.num_wpoints+s+num_blocks*blocksize+l]),
                                  snn[loffset]*n[i]*t[l]);
                    }
                }
            }
        }
    }
    return;
}


/*Calculate the optical depth contribution of the water vapor continuum.*/
__global__ void calc_water_vapor_ctm_optical_depth_d(uint64_t const num_wpoints,
                                                     int const num_layers,
                                                     fp_t * const tau,
                                                     fp_t const * const CS,
                                                     fp_t const * const T,
                                                     fp_t const * const Ps,
                                                     fp_t const * const N,
                                                     fp_t const * const T0,
                                                     fp_t const * const CF,
                                                     fp_t const * const P,
                                                     fp_t const * const T0F)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_wpoints)
    {
        fp_t const tref = 296.f;
        int i;
        for (i=0;i<num_layers;++i)
        {
            tau[i*num_wpoints+j] += N[i]*(tref/T[i])*((CS[j]*Ps[i]*
                                    EXP(T0[j]*(tref-T[i]))) +
                                    (CF[j]*(P[i]-Ps[i])*
                                    EXP(T0F[j]*(tref-T[i]))));
        }
    }
    return;
}


/*Calculate the optical depth contribution of the ozone continuum.*/
__global__ void calc_ozone_ctm_optical_depth_d(uint64_t const num_wpoints,
                                               int const num_layers,
                                               fp_t const * const cross_section,
                                               fp_t const * const N,
                                               fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_wpoints)
    {
        int i;
        for (i=0;i<num_layers;++i)
        {
            tau[i*num_wpoints+j] += N[i]*cross_section[j];
        }
    }
    return;
}


/*Do a quadratic interpolation of line wing values in each bin (except for the last one.).*/
__global__ void interpolate_d(SpectralBins_t const bins, fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < bins.n-1)
    {
        if (bins.do_interp)
        {
            int i;
            for (i=0;i<bins.num_layers;++i)
            {
                fp_t *t = &(tau[i*bins.num_wpoints]);
                fp_t const *x = &(bins.w[j*NIP]);
                fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
                bin_quad_interp(x,
                                y,
                                bins.l[j],
                                bins.r[j],
                                bins.w0,
                                bins.wres,
                                t);
            }
        }
        else
        {
            int i;
            for (i=0;i<bins.num_layers;++i)
            {
                fp_t *t = &(tau[i*bins.num_wpoints]);
                fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
                bin_no_interp(bins.l[j],
                              bins.r[j],
                              y,
                              t);
            }
        }
    }
    return;
}


/*Do a quadratic interpolation of line wing values in each bin (except for the last one.).*/
__global__ void interpolate_last_bin_d(SpectralBins_t const bins, fp_t * const tau)
{
    uint64_t const i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < bins.num_layers)
    {
        uint64_t const j = bins.n - 1;
        if (bins.do_last_interp)
        {
            fp_t *t = &(tau[i*bins.num_wpoints]);
            fp_t const *x = &(bins.w[j*NIP]);
            fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
            bin_quad_interp(x,
                            y,
                            bins.l[j],
                            bins.r[j],
                            bins.w0,
                            bins.wres,
                            t);
        }
        else
        {
            fp_t *t = &(tau[i*bins.num_wpoints]);
            fp_t const *y = &(bins.tau[i*bins.isize + j*NIP]);
            bin_no_interp(bins.l[j],
                          bins.r[j],
                          y,
                          t);
        }
    }
    return;
}


/*Calculate the optical depth contribution of a CFC.*/
__global__ void calc_cfc_optical_depth_d(uint64_t const num_wpoints, int const num_layers,
                                         fp_t const * const n, fp_t const * const x,
                                         fp_t const * const cross_section, fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_wpoints)
    {
        fp_t const half = 0.5;
        int i;
        for (i=0; i<num_layers; ++i)
        {
            tau[i*num_wpoints+j] += half*n[i]*(x[i]+x[i+1])*cross_section[j];
        }
    }
    return;
}


/*Calculate the optical depth contribution of collision-induced absorption.*/
__global__ void calc_cia_optical_depth_d(uint64_t const num_wpoints, int const num_layers,
                                         fp_t const * const p, fp_t const * const t,
                                         fp_t const * const x1, fp_t const * const x2,
                                         fp_t const * const cross_section, fp_t * const tau)
{
    uint64_t const j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j < num_wpoints)
    {
        fp_t const quarter = 0.25;
        fp_t const m = 28.97/6.02214076e23; /*[g].*/
        fp_t const g = 980.; /*[cm/s^2]*/
        fp_t const k = 1.38064852e-16; /*[(g*cm^2)/(s^2*K)].*/
        fp_t const atmtobarye = 1.013e6; /*[g/(cm*s^2*atm)].*/
        fp_t const c = (atmtobarye*atmtobarye)/(k*m*g*2.); /*[K/(atm^2*cm^5)].*/
        int i;
        for (i=0; i<num_layers; ++i)
        {
            fp_t n = c*((p[i]*p[i] - p[i+1]*p[i+1])/t[i])*quarter*(x1[i]+x1[i+1])*
                     (x2[i]+x2[i+1]);
            n = (n >= 0) ? n : n*-1.f;
            tau[i*num_wpoints+j] += n*cross_section[j];
        }
    }
    return;
}
#endif
