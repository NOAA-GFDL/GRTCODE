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

#include <stdio.h>
#include <stdlib.h>
#include "argparse.h"
#include "debug.h"
#include "netcdf.h"


#define nc_catch(e) { \
    int e_ = e; \
    if (e_ != NC_NOERR) {\
        fprintf(stderr, "[%s: %d] %s\n", __FILE__, __LINE__, nc_strerror(e_)); \
        exit(EXIT_FAILURE); \
    }}
#ifdef SINGLE_PRECISION
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_float(a, b, c, d, e))
#else
#define get_var(a, b, c, d, e) nc_catch(nc_get_vara_double(a, b, c, d, e))
#endif


int check_result(fp_t const * const actual, fp_t const * const expected, size_t const size,
                 fp_t const tolerance)
{
    int failures = 0;
    size_t i;
    for (i=0; i<size; ++i)
    {
        fp_t err = 100.*(fabs(actual[i] - expected[i]))/expected[i];
        if (err > tolerance)
        {
            failures++;
        }
    }
    return failures;
}


void read_data(int const ncid, char const * const varname, fp_t * const buffer,
               size_t const num_levels)
{
    int varid;
    nc_catch(nc_inq_varid(ncid, varname, &varid));
    size_t start[1] = {0};
    size_t count[1] = {num_levels};
    get_var(ncid, varid, start, count, buffer);
    return;
}


int main(int argc, char **argv)
{
    /*Add command line args.*/
    int const one = 1;
    Parser_t parser;
    parser = create_parser(argc, argv, NULL);
    add_argument(&parser, "input_file", NULL, "data file to check.", &one);
    add_argument(&parser, "-t", "--tolerance", "percent error tolerance.", &one);
    parse_args(parser);

    /*Open the file.*/
    char buf[valuelen];
    get_argument(parser, "input_file", buf);
    int ncid;
    nc_catch(nc_open(buf, NC_NOWRITE, &ncid));

    /*Get the number of levels and allocate a buffer.*/
    int dimid;
    nc_catch(nc_inq_dimid(ncid, "level", &dimid));
    size_t num_levels;
    nc_catch(nc_inq_dimlen(ncid, dimid, &num_levels));
    fp_t *buffer = NULL;
    gmalloc(buffer, num_levels, HOST_ONLY);

    /*Get the tolerance.*/
    fp_t tolerance = 1.;
    if (get_argument(parser, "-t", buf))
    {
        tolerance = atof(buf);
    }
    printf("Using a tolerance of %e%%\n", tolerance);
    int any_failures = 0;

    fp_t rlu[55] = {
        288.913737530876, 288.983737643993, 289.051925039948, 289.067916443805,
        288.987335448039, 288.846865788953, 288.719133322014, 288.585101779519,
        288.466404723276, 288.340862154713, 288.183852123831, 288.056668386618,
        288.025081922886, 288.032097059779, 287.988593791205, 287.914489313645,
        287.894254831224, 287.907515368633, 287.95840115505, 288.127259632008,
        288.270340559296, 288.477134219385, 289.062864717648, 289.930945225473,
        290.824066003184, 291.988567080749, 293.518002511827, 295.345653112413,
        297.527783502065, 299.953992830315, 302.573754494777, 305.246339024839,
        307.976144663881, 311.574533248324, 317.056265400995, 324.139873860774,
        331.699859241207, 338.669069571151, 345.361586521797, 352.133993834468,
        358.361457556774, 364.161006133891, 370.25898662835, 377.026588101515,
        384.363132495592, 391.492866335347, 397.805855689271, 402.942787915453,
        408.290952235571, 413.381956360222, 417.966934874191, 422.917460988538,
        428.39000272605, 433.120664868699, 437.254933154453};
    read_data(ncid, "rlu", buffer, num_levels);
    int failures = check_result(buffer, rlu, num_levels, tolerance);
    if (failures > 0)
    {
        printf("%d/%zu failures for upwelling longwave radiation.\n", failures, num_levels);
        any_failures += failures;
    }
    else
    {
        printf("All upwelling longwave radiation values pass.\n");
    }

    fp_t rld[55] = {
        0., 0.220737126841315, 0.392951503039547, 0.664264125049257,
        0.992828782327933, 1.30435309855765, 1.58288609606773, 1.88148170807822,
        2.18613475792047, 2.53616915363425, 2.92045577728199, 3.29421390241644,
        3.70438401479704, 4.26842988230056, 4.99651441280587, 5.78302195044713,
        6.5769154930942, 7.42191670576715, 8.2618540125893, 9.04911475530989,
        9.83063181766431, 10.334417016738, 10.4916187096612, 10.76610688554,
        11.3193980071403, 12.0533445321835, 13.2776386833051, 15.0996573196946,
        17.5410906148978, 20.8059472604357, 25.2273954555994, 31.3097078137587,
        41.3532425592645, 56.6122441374152, 73.8007349649491, 91.4694791696674,
        108.00769339831, 122.173154456698, 134.493398996982, 146.350471148304,
        157.753976202288, 169.183726371154, 182.182708702184, 196.652742370828,
        211.400674385671, 225.603168774383, 238.711724040232, 250.548340647059,
        261.445140262226, 272.662426479513, 283.042881984512, 292.486221801644,
        301.502641862324, 309.659907915117, 315.962577694042};
    read_data(ncid, "rld", buffer, num_levels);
    failures = check_result(buffer, rld, num_levels, tolerance);
    if (failures > 0)
    {
        printf("%d/%zu failures for downwelling longwave radiation.\n", failures, num_levels);
        any_failures += failures;
    }
    else
    {
        printf("All downwelling longwave radiation values pass.\n");
    }

    fp_t rsu[55] = {
        127.406821884485, 127.395538839041, 127.386035800316, 127.37453008599,
        127.361886398994, 127.349044301885, 127.337085455735, 127.327281425313,
        127.321101717635, 127.320129034226, 127.325881514769, 127.339700789144,
        127.362888782839, 127.396882811676, 127.442304449465, 127.497282203465,
        127.556996832682, 127.610574398576, 127.639035179513, 127.617434022047,
        127.518214408831, 127.312982850602, 126.988668610234, 126.539163792709,
        125.975760852428, 125.312718679304, 124.564164892749, 123.741999386037,
        122.856047727874, 121.913864554525, 120.921284857207, 119.883174698315,
        118.799968218914, 117.675436940302, 116.516587989312, 115.331256445066,
        114.123669059152, 112.89472421699, 111.649660469145, 110.397769800828,
        109.146006962174, 107.908759581337, 106.710505464018, 105.570394722533,
        104.494769646798, 103.485790255599, 102.551014808899, 101.705084291667,
        100.965154379498, 100.345931157551, 99.8535713702175, 99.4778925673346,
        99.1962679309566, 98.9887164955972, 98.8405181640125};
    read_data(ncid, "rsu", buffer, num_levels);
    failures = check_result(buffer, rsu, num_levels, tolerance);
    if (failures > 0)
    {
        printf("%d/%zu failures for upwelling shortwave radiation.\n", failures, num_levels);
        any_failures += failures;
    }
    else
    {
        printf("All upwelling shortwave radiation values pass.\n");
    }

    fp_t rsd[55] = {
        757.354656462478, 757.103209175312, 756.800438451638, 756.251500253976,
        755.456311827089, 754.473376920881, 753.393540776989, 752.30185747332,
        751.237272125942, 750.184409517985, 749.10609543696, 747.971624099202,
        746.765747816562, 745.478555518823, 744.099310735927, 742.62431338666,
        741.053813749215, 739.406582680514, 737.720536046163, 736.034660525226,
        734.374483575057, 732.748678196562, 731.238074187011, 729.815334476592,
        728.426247337693, 727.003149694526, 725.509667438775, 723.933766686086,
        722.260025978381, 720.465244739009, 718.453744457607, 716.069265583608,
        712.496740524299, 706.335862093769, 697.896014622244, 688.606955886987,
        679.913634556508, 672.557891408052, 666.176915486093, 660.296622973967,
        654.92643179964, 649.543162717245, 643.277484849939, 636.159046381035,
        628.907942429397, 622.017008524741, 615.528752744182, 609.285934814703,
        603.156976269595, 597.129325040076, 591.323208415101, 585.982928882095,
        581.324549805146, 577.40324657116, 574.177651956291};
    read_data(ncid, "rsd", buffer, num_levels);
    failures = check_result(buffer, rsd, num_levels, tolerance);
    if (failures > 0)
    {
        printf("%d/%zu failures for downwelling shortwave radiation.\n", failures, num_levels);
        any_failures += failures;
    }
    else
    {
        printf("All downwelling shortwave radiation values pass.\n");
    }
    nc_catch(nc_close(ncid));
    gfree(buffer, HOST_ONLY);
    destroy_parser(&parser);
    return any_failures;
}
