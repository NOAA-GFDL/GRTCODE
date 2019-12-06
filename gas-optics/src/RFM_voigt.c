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


/* This is the voishp subroutine from the Reference Forward Model written by A. Dudhia
   and ported to c by Raymond Menzel.  A paper describing the Reference Forward Model
   is located at https://doi.org/10.1016/j.jqsrt.2016.06.018.

  VERSION
      19-SEP-11  AD  Move DATA statements after variable declarations
      03-MAY-00  RJW JQSRT corrections
      30-DEC-98  RJW Incorporate JQSRT paper version of Humlicek revision
      03-MAR-97  AD  Version 3.
      01-OCT-96  AD  Version 2.
      01-SEP-96  AD  Version 1.
      16-JUL-96  RJW Corrected typo
      02-JUL-96  RJW Modified Humlicek region 1 (assumed major region)
      08-JUN-96  AD  Original. Based on GENLN2 module VOIGT.

  DESCRIPTION
      Calculate Voigt Line shape.
      Called by RFMFIN and RFMWID.
      This modifies the GENLN2 VOIGT algorithm - see VOIGL2 for original.
      The Voigt lineshape formulation:

                  g(X,Y) = S * g0 * K(X,Y)
                  g0 = 1/Ad * SQRT(ln2/pi)
                  X = (nu - nu0)/Ad *SQRT(ln2)
                  Y = Al/Ad *SQRT(ln2)
                  K(X,Y) = Y/pi *
                  INT^(+infty)_(-infty){exp(-t**2)/[Y**2 + (X-t)**2]}dt

      This routine calculates the complex probability function using a
      modified version of the Humlicek algorithm (JQSRT V27 437 1982)
      accepted for publication in JQSRT 1999.
      The calculation is performed for the array of x,y pairs for a given line
      over the fine mesh points of the current wide mesh.
*/


#include <math.h>
#include <stdint.h>
#include "debug.h"
#include "grtcode_config.h"
#include "grtcode_utilities.h"
#include "line_shape.h"
#include "RFM_voigt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/*pi^(-1/2).*/
#ifdef RSQRPI
#error
#else
#define RSQRPI 0.56418958f
#endif

/*(ln(2))^(1/2).*/
#ifdef SQRLN2
#error
#else
#define SQRLN2 0.832554611f
#endif


/*Calculate the Voigt line shape function using the Humlicek algorithm, as
  implemented in the Reference Forward Model (RFM).*/
HOST DEVICE int rfm_voigt_line_shape(LineShapeInputs_t const vals, fp_t * const K)
{
    fp_t const DWNO = vals.w;
    fp_t const WNOADJ = vals.line_center;
    fp_t const WIDADJ = vals.lorentz_hwhm;
    fp_t const DOPADJ = vals.doppler_hwhm;
    fp_t const WRES = vals.wres;
    uint64_t const N = vals.num_wpoints;

    float const REPWID = SQRLN2/DOPADJ;
    float const Y = REPWID*WIDADJ;
    float const YQ = Y*Y;
    if (Y >= 70.55f)
    {
        uint64_t i;
        for (i=0;i<N;++i)
        {
            float const XI = (DWNO + i*WRES - WNOADJ)*REPWID;
            K[i] = REPWID*Y/(M_PI*(XI*XI+YQ));
        }
        return GRTCODE_SUCCESS;
    }

    float const YRRTPI = Y*RSQRPI;
    float const XLIM0 = SQRT(15100.0f + Y*(40.0f - Y*3.6f));
    float XLIM1;
    if (Y >= 8.425f)
    {
        XLIM1 = 0.0f;
    }
    else
    {
        XLIM1 = SQRT(164.0f - Y*(4.3f + Y*1.8f));
    }
    float XLIM2 = 6.8f - Y;
    float const XLIM3 = 2.4f*Y;
    float const XLIM4 = 18.1f*Y + 1.65f;
    if (Y <= 0.000001f)
    {
        XLIM1 = XLIM0;
        XLIM2 = XLIM0;
    }

    int RG1 = 1;
    float A0;
    float D0;
    float D2;
    int RG2 = 1;
    float H0;
    float H2;
    float H4;
    float H6;
    float E0;
    float E2;
    float E4;
    int RG3 = 1;
    float Z0;
    float Z2;
    float Z4;
    float Z6;
    float Z8;
    float P0;
    float P2;
    float P4;
    float P6;
    float P8;
    float const Y0 = 1.5f;
    float const Y0PY0 = 3.f; 
    float const Y0Q = 2.25f;
    float const YPY0 = Y + Y0;
    float const YPY0Q = YPY0*YPY0;
    float const C[6] = {1.0117281f, -0.75197147f, 0.012557727f,
                        0.010022008f, -0.00024206814f, 0.00000050084806f};
    float const S[6] = {1.393237f, 0.23115241f, -0.15535147f,
                        0.0062183662f, 0.000091908299f, -0.00000062752596f};
    float const T[6] = {0.31424038f, 0.94778839f, 1.5976826f,
                        2.2795071f, 3.0206370f, 3.8897249f};
    uint64_t i;
    for (i=0;i<N;++i)
    {
        float const XI = (DWNO + i*WRES - WNOADJ)*REPWID;
        float const ABX = ABS(XI);
        float const XQ = ABX*ABX;
        if (ABX >= XLIM0)
        {
            K[i] = YRRTPI/(XQ + YQ);
        }
        else if (ABX >= XLIM1)
        {
            if (RG1 != 0)
            {
                RG1 = 0;
                A0 = YQ + 0.5;
                D0 = A0*A0;
                D2 = YQ + YQ - 1.0;
            }
            float const D = RSQRPI/(D0 + XQ*(D2 + XQ));
            K[i] = D*Y*(A0 + XQ);
        }
        else if (ABX >= XLIM2)
        {
            if (RG2 != 0)
            {
                RG2 = 0;
                H0 = 0.5625f + YQ*(4.5f + YQ*(10.5f + YQ*(6.0f + YQ)));
                H2 = -4.5f + YQ*(9.0f + YQ*(6.0f + YQ*4.0f));
                H4 = 10.5f - YQ*(6.0f - YQ*6.0f);
                H6 = -6.0f + YQ* 4.0f;
                E0 = 1.875f + YQ*(8.25f + YQ*(5.5f + YQ));
                E2 = 5.25f + YQ*(1.0f + YQ*3.0f);
                E4 = 0.75f*H6;
            }
            float const D = RSQRPI/(H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))));
            K[i] = D*Y*(E0 + XQ*(E2 + XQ*(E4 + XQ)));
        }
        else if (ABX < XLIM3)
        {
            if (RG3 != 0)
            {
                RG3 = 0;
                Z0 = 272.1014f + Y*(1280.829f + Y*(2802.870f + Y*(3764.966f
                     + Y*(3447.629f + Y*(2256.981f + Y*(1074.409f + Y*(369.1989f
                     + Y*(88.26741f + Y*(13.39880f + Y)))))))));
                Z2 = 211.678f + Y*(902.3066f + Y*(1758.336f + Y*(2037.310f
                     + Y*(1549.675f + Y*(793.4273f + Y*(266.2987f
                     + Y*(53.59518f + Y*5.0f)))))));
                Z4 = 78.86585f + Y*(308.1852f + Y*(497.3014f + Y*(479.2576f
                     + Y*(269.2916f + Y*(80.39278f + Y*10.0f)))));
                Z6 = 22.03523f + Y*(55.02933f + Y*(92.75679f + Y*(53.59518f
                     + Y*10.0f)));
                Z8 = 1.496460f + Y*(13.39880f + Y*5.0f);
                P0 = 153.5168f + Y*(549.3954f + Y*(919.4955f + Y*(946.8970f
                     + Y*(662.8097f + Y*(328.2151f + Y*(115.3772f + Y*(27.93941f
                     + Y*(4.264678f + Y*0.3183291f))))))));
                P2 = -34.16955f + Y*(-1.322256f + Y*(124.5975f + Y*(189.7730f
                     + Y*(139.4665f + Y*(56.81652f + Y*(12.79458f
                     + Y*1.2733163f))))));
                P4 = 2.584042f + Y*(10.46332f + Y*(24.01655f + Y*(29.81482f
                     + Y*(12.79568f + Y*1.9099744f))));
                P6 = -0.07272979f + Y*(0.9377051f + Y*(4.266322f + Y*1.273316f));
                P8 = 0.0005480304f + Y*0.3183291f;
            }
            float const D = 1.7724538f/(Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 +
                            XQ*(Z8+XQ)))));
            K[i] = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))));
        }
        else
        {
            K[i] = 0.0f;
            float MQ[6];
            float MF[6];
            float XM[6];
            float YM[6];
            float PQ[6];
            float PF[6];
            float XP[6];
            float YP[6];
            int J;
            for (J=0;J<=5;J++)
            {
                float D = XI - T[J];
                MQ[J] = D*D;
                MF[J] = 1.0f/(MQ[J] + YPY0Q);
                XM[J] = MF[J]*D;
                YM[J] = MF[J]*YPY0;
                D = XI + T[J];
                PQ[J] = D*D;
                PF[J] = 1.0f/(PQ[J] + YPY0Q);
                XP[J] = PF[J]*D;
                YP[J] = PF[J]*YPY0;
            }

            if (ABX <= XLIM4)
            {
                for (J=0;J<=5;J++)
                {
                    K[i] = K[i] + C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]);
                }
            }
            else
            {
                float const YF = Y + Y0PY0;
                for (J=0;J<=5;J++)
                {
                    K[i] = K[i]
                           + (C[J]*(MQ[J]*MF[J]-Y0*YM[J]) + S[J]*YF*XM[J])/
                           (MQ[J]+Y0Q)
                           + (C[J]*(PQ[J]*PF[J]-Y0*YP[J]) - S[J]*YF*XP[J])/
                           (PQ[J]+Y0Q);
                }
                K[i] = Y*K[i] + EXP(-XQ);
            }
        }
        K[i] = RSQRPI*REPWID*K[i];
    }
    return GRTCODE_SUCCESS;
}
