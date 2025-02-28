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

#ifndef DEBUG_H_
#define DEBUG_H_

#include <errno.h>
#include <fenv.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "device.h"
#include "return_codes.h"
#include "utilities.h"
#include "verbosity.h"
#include "verbosity-internal.h"


#define backtrace() { \
    char s____[1024]; \
    snprintf(s____, 1024, "\r\33[2K\t%s: %d\n", __FILE__, __LINE__); \
    append_to_error_buffer(s____); \
}


#define log_warn(mesg, ...) { \
    if (grtcode_verbosity() >= GRTCODE_WARN) { \
        char s__[1024]; \
        snprintf(s__, 1024, mesg, __VA_ARGS__); \
        fprintf(stderr, "\r\33[2K[%s:%d] warning: %s\n", __FILE__, __LINE__, s__); \
    } \
}


#define log_info(mesg, ...) { \
    if (grtcode_verbosity() >= GRTCODE_INFO) { \
        char s__[1024]; \
        snprintf(s__, 1024, mesg, __VA_ARGS__); \
        fprintf(stderr, "\r\33[2K[%s:%d] info: %s\n", __FILE__, __LINE__, s__); \
    } \
}


#define log_mesg(mesg, ...) { \
    if (grtcode_verbosity() >= GRTCODE_NONE) { \
        char s__[1024]; \
        snprintf(s__, 1024, mesg, __VA_ARGS__); \
        fprintf(stdout, "\r\33[2K %s\n", s__); \
    } \
}


/*Macros that return error codes.*/
#define raise(err, mesg, ...) { \
    reset_error_buffer(); \
    char s__[1024]; \
    snprintf(s__, 1024, mesg , __VA_ARGS__); \
    char s___[2048]; \
    snprintf(s___, 2048, "Error: %s\nBacktrace:", s__); \
    append_to_error_buffer(s___); \
    backtrace(); \
    return err; \
}


#ifdef __CUDA_ARCH__
#define catch(val) { \
    int e_ = val; \
    if (e_ != GRTCODE_SUCCESS) { \
        return e_; \
    } \
}
#else
#define catch(val) { \
    int e_ = val; \
    if (e_ != GRTCODE_SUCCESS) { \
        backtrace(); \
        return e_; \
    } \
}
#endif


/*Safety checks.*/
#if defined(__CUDA_ARCH__) || defined(FAST)
#define sentinel() {return GRTCODE_SENTINEL_ERR;}
#define not_null(p) {}
#define is_null(p) {}
#define not_nan(v) {}
#define min_check(v, min) {}
#define max_check(v, max) {}
#define in_range(v, min, max) {}
#define assert(v1, v2) {}
#define clear_floating_point_exceptions() {}
#define floating_point_error_code(e) {}
#define catch_floating_point_exceptions(e) {}
#else
#define sentinel() { \
    char const *s_ = "This branch should never be reached (%s,%d)."; \
    raise(GRTCODE_SENTINEL_ERR, s_, __FILE__, __LINE__); \
}


#define not_null(p) { \
    if (p == NULL) { \
        char const *s_ = "null pointer at address %p."; \
        raise(GRTCODE_NULL_ERR, s_, (void *)(&p)); \
    } \
}


#define is_null(p) { \
    if (p != NULL) { \
        char const *s_ = "pointer at address %p is not null."; \
        raise(GRTCODE_NON_NULL_ERR, s_, (void *)(&p)); \
    } \
}


#define not_nan(v) { \
    if (isnan((double)v)) { \
        char const *s_ = "input value (%e) is Nan."; \
        raise(GRTCODE_INVALID_ERR, s_, (double)v); \
    } \
}


#define min_check(v, min) { \
    not_nan(v); \
    not_nan(min); \
    if (v < min) { \
        char const *s_ = "value (%e) less than minimum allowed (%e)."; \
        raise(GRTCODE_RANGE_ERR, s_, (double)v, (double)min); \
    } \
}


#define max_check(v, max) { \
    not_nan(v); \
    not_nan(max); \
    if (v > max) { \
        char const *s_ = "value (%e) greater than maximum allowed (%e)."; \
        raise(GRTCODE_RANGE_ERR, s_, (double)v, (double)max); \
    } \
}


#define in_range(v, min, max) { \
    min_check(v, min); \
    max_check(v, max); \
    if (min > max) { \
        char const *s_ = "min value (%e) greater tha max value (%e)."; \
        raise(GRTCODE_RANGE_ERR, s_, (double)min, (double)max); \
    } \
}


#define assert(v1, v2) { \
    if (v1 != v2) { \
        char const *s_ = "values (%u, %u) are not equal."; \
        raise(GRTCODE_VALUE_ERR, s_, (uint64_t)v1, (uint64_t)v2); \
    } \
}


#if defined(FE_DIVBYZERO) && defined(FE_INEXACT) && defined(FE_INVALID) \
    && defined(FE_OVERFLOW) && defined(FE_UNDERFLOW)
/*
#pragma STDC FENV_ACCESS ON
*/


#define clear_floating_point_exceptions() { \
    if (math_errhandling & MATH_ERREXCEPT) { \
        feclearexcept(FE_ALL_EXCEPT); \
    } \
    errno = 0; \
}


#define floating_point_error_code(e) { \
    e = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW); \
    if (e == 0) { \
        e = GRTCODE_SUCCESS; \
    } \
    else if ((math_errhandling & MATH_ERREXCEPT) && e) { \
        if (FE_DIVBYZERO & e) { \
            e = GRTCODE_DIVBYZERO_ERR; \
        } \
        else if (FE_OVERFLOW & e) { \
            e = GRTCODE_OVERFLOW_ERR; \
        } \
        else if (FE_INVALID & e) { \
            e = GRTCODE_INVALID_ERR; \
        } \
        else if (FE_UNDERFLOW & e) { \
            e = GRTCODE_SUCCESS; \
        } \
        else { \
            sentinel(); \
        } \
    } \
    else if ((math_errhandling & MATH_ERRNO) && errno != 0) \
    { \
        e = GRTCODE_INVALID_ERR; \
    } \
    else { \
        sentinel(); \
    } \
}


#define catch_floating_point_exceptions() { \
    int e_; \
    floating_point_error_code(e_); \
    if (e_ == GRTCODE_DIVBYZERO_ERR) { \
        char const *mesg = "divide by zero (fetestexcept=%d)."; \
        raise(e_, mesg, FE_DIVBYZERO); \
    } \
    else if (e_ == GRTCODE_OVERFLOW_ERR) { \
        char const *mesg = "overflow (fetestexcept=%d)."; \
        raise(e_, mesg, FE_OVERFLOW); \
    } \
    else if (e_ == GRTCODE_INVALID_ERR) { \
        char const *mesg = "floating point invalid (fetestexcept=%d)."; \
        raise(e_, mesg, FE_INVALID); \
    } \
}
#else
#pragma message ("Floating point exception handling is not used.")
#define clear_floating_point_exceptions() {}
#define floating_point_error_code(e) {}
#define catch_floating_point_exceptions(e) {}
#endif


#endif


#define cat(a,b) a##b
#define str(a) #a


#ifdef __NVCC__
#define gpu_catch(val) { \
    cudaError_t e_ = val; \
    if (e_ != cudaSuccess) { \
        char const *s_ = "cuda: %s"; \
        raise(GRTCODE_GPU_ERR, s_, cudaGetErrorString(e_)); \
    } \
}


#define _glaunch(func, threads, loc,...) { \
    gpu_catch(cudaSetDevice(loc)); \
    int min_grid_size; \
    int dim_block; \
    gpu_catch(cudaOccupancyMaxPotentialBlockSize(&min_grid_size, &dim_block, cat(func, _d), \
                                                 0, (int)threads)); \
    int dim_grid = (((int)threads) + dim_block - 1)/dim_block; \
    log_info("Running %s with %d thread-blocks and %d threads-per-thread-block on device %d.", \
             str(func), dim_grid, dim_block, loc); \
    cat(func, _d)<<<dim_grid, dim_block, 0, 0>>>(__VA_ARGS__); \
}


#define FROM_HOST cudaMemcpyHostToDevice
#define FROM_DEVICE cudaMemcpyDeviceToHost
#define HOST __host__
#define DEVICE __device__
#else
#define gpu_catch(val) {}
#define _glaunch(func, threads, loc, ...) {}
#define FROM_HOST
#define FROM_DEVICE
#define HOST
#define DEVICE
#endif


#ifndef _OPENMP
#define omp_set_num_threads(n) {}
#define omp_get_max_threads() 1
#endif


#define gmalloc(ptr, size, loc) { \
    if (loc == HOST_ONLY) { \
        catch(malloc_ptr((void **)&ptr, sizeof(*ptr)*size)); \
    } \
    else { \
        gpu_catch(cudaSetDevice(loc)); \
        gpu_catch(cudaMalloc(&ptr, sizeof(*ptr)*size)); \
    } \
}


#define gfree(ptr, loc) { \
    if (loc == HOST_ONLY) { \
        catch(free_ptr((void **)&ptr)); \
    } \
    else { \
        gpu_catch(cudaSetDevice(loc)); \
        gpu_catch(cudaFree(ptr)); \
    } \
}


#define gmemset(ptr, val, size, loc) { \
    if (loc == HOST_ONLY) { \
        memset(ptr, val, sizeof(*ptr)*size); \
    } \
    else { \
        gpu_catch(cudaSetDevice(loc)); \
        gpu_catch(cudaMemset(ptr, val, sizeof(*ptr)*size)); \
    } \
}


#define gmemcpy(dst, src, size, loc, dir) { \
    if (loc == HOST_ONLY) { \
        memcpy(dst, src, sizeof(*dst)*size); \
    } \
    else { \
        gpu_catch(cudaSetDevice(loc)); \
        gpu_catch(cudaMemcpy(dst, src, sizeof(*dst)*size, dir)); \
    } \
}


#define glaunch(func, threads, loc, ...) { \
    if ((int)threads > 0) \
    { \
        if (loc == HOST_ONLY) { \
            int t_ = omp_get_max_threads(); \
            if ((uint64_t)threads < (uint64_t)t_) { \
                omp_set_num_threads(threads); \
            } \
            log_info("Running %s with %d openmp threads on host.", str(func), \
                     ((uint64_t)threads < (uint64_t)t_) ? (int)threads : t_); \
            catch(func(__VA_ARGS__)); \
            if ((uint64_t)threads < (uint64_t)t_) { \
                omp_set_num_threads(t_); \
            } \
        } \
        else { \
            _glaunch(func, threads, loc, __VA_ARGS__); \
        } \
    } \
}


#endif
