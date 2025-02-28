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

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "parse_HITRAN_file.h"
#include "parse_HITRAN_file-internal.h"
#include "tips2017.h"


typedef enum HITRAN2012_cols
{
    mol_c,
    iso_c,
    Vnn_c,
    Snn_c,
    A_c,
    Yair_c,
    Yself_c,
    Elo_c,
    n_c,
    del_c,
    Vu_c,
    Vl_c,
    Qu_c,
    Ql_c,
    Ierr_c,
    Iref_c,
    flag_c,
    gu_c,
    gl_c,
    NCOLS
} HITRAN2012_col_t;


typedef enum LookupCast
{
    NIL,
    I32,
    F32,
    F64
} LookupCast_t;


typedef union HITRAN2012_vals
{
    void *nil;
    int i;
    float f;
    double d;
} HITRAN2012_vals_t;


static unsigned int const HITRAN2012_recordLen = 160;
static unsigned int const HITRAN2012_pad = 2;
static unsigned int const HITRAN2012_fmt[NCOLS][2] =
{
    {2 , I32},
    {1 , I32},
    {12, F64},
    {10, F64},
    {10, NIL},
    {5 , F32},
    {5,  F32},
    {10, F32},
    {4,  F32},
    {8,  F32},
    {15, NIL},
    {15, NIL},
    {15, NIL},
    {15, NIL},
    {6,  NIL},
    {12, NIL},
    {1,  NIL},
    {7,  NIL},
    {7,  NIL}
};


/*Allocate memory for molecular line parameters.*/
static int alloc_line_params(LineParams_t * const line_params, /*Molecular line parameters.*/
                             uint64_t const num_lines, /*Number of molecular lines.*/
                             Device_t const device /*Device id.*/
                            )
{
    not_null(line_params);
    line_params->num_lines = num_lines;
    line_params->device = device;
    gmalloc(line_params->iso, num_lines, device);
    gmalloc(line_params->vnn, num_lines, device);
    gmalloc(line_params->snn, num_lines, device);
    gmalloc(line_params->yair, num_lines, device);
    gmalloc(line_params->yself, num_lines, device);
    gmalloc(line_params->en, num_lines, device);
    gmalloc(line_params->n, num_lines, device);
    gmalloc(line_params->d, num_lines, device);
    return GRTCODE_SUCCESS;
}


/*Free memory used for molecular line parameters.*/
int free_line_params(LineParams_t * const line_params)
{
    not_null(line_params);
    gfree(line_params->iso, line_params->device);
    gfree(line_params->vnn, line_params->device);
    gfree(line_params->snn, line_params->device);
    gfree(line_params->yair, line_params->device);
    gfree(line_params->yself, line_params->device);
    gfree(line_params->en, line_params->device);
    gfree(line_params->n, line_params->device);
    gfree(line_params->d, line_params->device);
    return GRTCODE_SUCCESS;
}


/*Reallocate the use only a subset of all lines.*/
static int realloc_line_params(LineParams_t * const line_params, /*Molecular line parameters.*/
                               uint64_t const num_lines /*Number of subset of lines.*/
                              )
{
    not_null(line_params);
    LineParams_t t;
    gmemcpy(&t, line_params, 1, HOST_ONLY, FROM_HOST);
    catch(alloc_line_params(line_params, num_lines, t.device));
    gmemcpy(line_params->iso, t.iso, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->vnn, t.vnn, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->snn, t.snn, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->yair, t.yair, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->yself, t.yself, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->en, t.en, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->n, t.n, num_lines, t.device, FROM_HOST);
    gmemcpy(line_params->d, t.d, num_lines, t.device, FROM_HOST);
    catch(free_line_params(&t));
    return GRTCODE_SUCCESS;
}


/*Cast the value from the HITRAN input file to the correct type.*/
static int HITRAN2012_cast(HITRAN2012_vals_t * const val, /*Output value.*/
                           int const col, /*Column index.*/
                           char const * const sval /*Input value.*/
                          )
{
    not_null(val);
    not_null(sval);
    LookupCast_t const typ = (LookupCast_t)(HITRAN2012_fmt[col][1]);
    switch (typ)
    {
        case NIL:
            val->nil = NULL;
            break;
        case I32:
            if (col == 1)
            {
                /*Handle HITRAN edge cases.  Since only 1 character is used
                  to specify the isotopologue, but the number of
                  isotopologues can exceed 9, they use assume the character
                  is hex-like (e.g., 0 = 10, A = 11, B = 12, ...).*/
                char c = sval[0];
                if (c == '0')
                {
                    val->i = 10;
                    break;
                }
                if (c >= 'A' && c <= 'Z')
                {
                    val->i = c - 'A' + 11;
                    break;
                }
            }
            catch(to_int(sval, &(val->i)));
            break;
        case F64:
        case F32:
            catch(to_double(sval, &(val->d)));
            if (typ == F32)
            {
                if (val->d >= -1.f*FLT_MAX && val->d <= FLT_MAX)
                {
                    val->f = val->d;
                }
                else
                {
                    char const *mesg = "value %e from column %d cannot be safely"
                                       " cast as a float.";
                    raise(GRTCODE_VALUE_ERR, mesg, val->d, col);
                }
            }
            break;
        default:
            {char const *mesg = "cast failed on col %d, LookupCast_t %d, sval: %s.";
            raise(GRTCODE_VALUE_ERR, mesg, HITRAN2012_fmt[col][0], HITRAN2012_fmt[col][1],
                  sval);}
    }
    return GRTCODE_SUCCESS;
}


/*Read molecular line parameters from HITRAN file.*/
int parse_hitran_file(LineParams_t * const line_params, char const * const filepath,
                      int const mol_id, double const w0, double const wn, Device_t const device)
{
    not_null(line_params);
    not_null(filepath);

    /*Open the file.*/
    char const *mesg = "Opening and reading HITRAN line parameters from file %s.";
    log_info(mesg, filepath);
    FILE *fp = NULL;
    catch(open_file(&fp, filepath, "r"));

    /*Count the number of lines in the file.*/
    size_t const max_line = 163;
    char *buf;
    gmalloc(buf, max_line, HOST_ONLY);
    gmemset(buf, 0, max_line, HOST_ONLY);
    ssize_t ll = 0;
    size_t l = 0;
    uint64_t n = 0;
    while ((ll = getline(&buf, &l, fp)) != -1)
    {
        ++n;
    }
    rewind(fp);

    /*Malloc space.*/
    LineParams_t lp;
    catch(alloc_line_params(&lp, n, HOST_ONLY));

    /*Parse out the line parameters.*/
    n = 0;
    size_t line_count = 0;
    while((ll = getline(&buf, &l, fp)) != -1)
    {
        line_count++;
        if ((ll - HITRAN2012_recordLen) > HITRAN2012_pad)
        {
            mesg = "Found bad record at line %zu (%zu exceeds max %zu chars)"
                   " in file %s.";
            raise(GRTCODE_VALUE_ERR, mesg, line_count, (size_t)ll, max_line, filepath);
        }
        else if (ll < HITRAN2012_recordLen)
        {
            mesg = "Found bad record at line %zu (%zu less than %zu chars)"
                   " in file %s.";
            raise(GRTCODE_VALUE_ERR, mesg, line_count, (size_t)ll, max_line, filepath);
        }
        unsigned int val_idx = 0;
        size_t offset = 0;
        int col;
        int go_to_next_line = 0;
        for (col=0; col<NCOLS; ++col)
        {
            if (go_to_next_line)
            {
                break;
            }
            size_t len = HITRAN2012_fmt[col][0];
            unsigned int t = HITRAN2012_fmt[col][1];
            char tmp[16];
            strncpy(tmp, &(buf[offset]), len);
            tmp[len]='\0';
            offset += len;
            if (t != NIL)
            {
                HITRAN2012_vals_t val;
                catch(HITRAN2012_cast(&val, col, tmp));

                /*Column indices for parsing.*/
                enum RefLinePtrIdx
                {
                    mol_pidx,
                    iso_pidx,
                    Vnn_pidx,
                    Snn_ref_pidx,
                    Yair_pidx,
                    Yself_pidx,
                    En_pidx,
                    n_pidx,
                    d_pidx
                };

                switch (val_idx)
                {
                    case mol_pidx:
                        if (val.i != mol_id)
                        {
                            go_to_next_line = 1;
                            continue;
                        }
                        break;
                    case iso_pidx:
                        lp.iso[n] = val.i;
                        break;
                    case Vnn_pidx:
                        lp.vnn[n] = (fp_t)(val.d);
                        break;
                    case Snn_ref_pidx:
                        lp.snn[n] = (fp_t)(val.d);
                        break;
                    case Yair_pidx:
                        lp.yair[n] = val.f;
                        break;
                    case Yself_pidx:
                        lp.yself[n] = val.f;
                        break;
                    case En_pidx:
                        lp.en[n] = val.f;
                        break;
                    case n_pidx:
                        lp.n[n] = val.f;
                        break;
                    case d_pidx:
                        lp.d[n] = val.f;
                        break;
                    default:
                        mesg = "Unknown column index (%d) on line %zu in file %s.";
                        raise(GRTCODE_VALUE_ERR, mesg, val_idx, line_count, filepath);
                }
                ++val_idx;
            }
        }
        if (!go_to_next_line)
        {
            if ((w0 < 0 && wn < 0) || (lp.vnn[n] >= w0 && lp.vnn[n] <= wn))
            {
                /*Include the line for the calculation.*/
                ++n;
            }
        }
    }
    gfree(buf, HOST_ONLY);

    /*Close the file.*/
    if (fclose(fp))
    {
        mesg = "error closing file %s.";
        raise(GRTCODE_IO_ERR, mesg, filepath);
    }

    /*Reallocate if necessary.*/
    if (lp.num_lines != n)
    {
        catch(realloc_line_params(&lp, n));
    }

    /*Adjust the raw read-in line strengths.*/
    fp_t const tref = 296.f;
    fp_t const c2 = -1.4387686f;
    fp_t *snn = lp.snn;
    int const *iso = lp.iso;
    fp_t const *en = lp.en;
    fp_t const *vnn = lp.vnn;
    uint64_t i;
#pragma omp parallel for default(shared) private(i)
    for (i=0; i<n; ++i)
    {
        snn[i] *= Q(mol_id, tref, iso[i])/(EXP(c2*en[i]/tref)*
                  (1.f - EXP(c2*vnn[i]/tref)));
    }

    if (device == HOST_ONLY)
    {
        line_params->num_lines = lp.num_lines;
        line_params->device = lp.device;
        line_params->iso = lp.iso;
        line_params->vnn = lp.vnn;
        line_params->snn = lp.snn;
        line_params->yair = lp.yair;
        line_params->yself = lp.yself;
        line_params->en = lp.en;
        line_params->n = lp.n;
        line_params->d = lp.d;
    }
    else
    {
        catch(alloc_line_params(line_params, lp.num_lines, device));
        gmemcpy(line_params->iso, lp.iso, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->vnn, lp.vnn, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->snn, lp.snn, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->yair, lp.yair, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->yself, lp.yself, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->en, lp.en, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->n, lp.n, lp.num_lines, device, FROM_HOST);
        gmemcpy(line_params->d, lp.d, lp.num_lines, device, FROM_HOST);
        catch(free_line_params(&lp));
    }
    return GRTCODE_SUCCESS;
}
