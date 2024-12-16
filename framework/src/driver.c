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
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "argparse.h"
#include "cloud_optics.h"
#include "distribute.h"
#include "driver.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "longwave.h"
#include "rayleigh.h"
#include "shortwave.h"
#include "solar_flux.h"

#ifdef _OPENMP
#include "omp.h"
#endif


#define catch(e) { \
    int e_ = e; \
    if (e_ != GRTCODE_SUCCESS) { \
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        char b_[1024]; \
        grtcode_errstr(e_, b_, 1024); \
        fprintf(stderr, "%s", b_); \
        return EXIT_FAILURE; \
    }}


/*Integrate using simple trapezoids on a uniform grid.*/
static void integrate(fp_t const * const in, /*Data to be integrated.*/
                      uint64_t const n, /*Size of input data array.*/
                      fp_t const dx, /*Grid spacing.*/
                      fp_t * const out /*Result of integral.*/
                     )
{
    *out = 0.;
    uint64_t i;
    for (i=0; i<n-1; ++i)
    {
        *out += dx*0.5*(in[i] + in[i+1]);
    }
    return;
}


/*Add molecules, CFCs, and collision-induced absoprtion.*/
static int add_molecules(GasOptics_t * const lbl, /*Gas optics object.*/
                         Atmosphere_t const atm /*Atmospheric state.*/
                        )
{
    int i;
    for (i=0; i<atm.num_molecules; ++i)
    {
        catch(add_molecule(lbl, atm.molecules[i], NULL, NULL));
    }
    for (i=0; i<atm.num_cfcs; ++i)
    {
        catch(add_cfc(lbl, atm.cfc[i].id, atm.cfc[i].path));
    }
    for (i=0; i<atm.num_cias; ++i)
    {
        catch(add_cia(lbl, atm.cia[i].id[0], atm.cia[i].id[1], atm.cia[i].path));
    }
    return GRTCODE_SUCCESS;
}


/*Calculate aerosol optics.*/
static int calculate_aerosol_optics(Atmosphere_t const atm, /*Atmospheric state.*/
                                    int const column, /*Column index.*/
                                    int const offset, /*Additional index offset.*/
                                    Optics_t * const aerosol /*Optics due to aerosols.*/
                                   )
{
    int i = column + offset;
    int n = atm.num_layers*aerosol->grid.n;
    fp_t *tau = (fp_t *)malloc(sizeof(*tau)*n);
    fp_t *omega = (fp_t *)malloc(sizeof(*omega)*n);
    fp_t *g = (fp_t *)malloc(sizeof(*g)*n);
    int j;
    for (j=0; j<atm.num_layers; ++j)
    {
        n = i*atm.num_layers*atm.aerosol_grid_size + j*atm.aerosol_grid_size;
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_optical_depth[n]),
                                  atm.aerosol_grid_size, &(tau[j*aerosol->grid.n]),
                                  linear_sample, NULL));
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_single_scatter_albedo[n]),
                                  atm.aerosol_grid_size, &(omega[j*aerosol->grid.n]),
                                  linear_sample, NULL));
        catch(interpolate_to_grid(aerosol->grid, atm.aerosol_grid, &(atm.aerosol_asymmetry_factor[n]),
                                  atm.aerosol_grid_size, &(g[j*aerosol->grid.n]),
                                  linear_sample, NULL));
    }
    catch(update_optics(aerosol, tau, omega, g));
    free(tau);
    free(omega);
    free(g);
    return GRTCODE_SUCCESS;
}


/*Calculate gas optics.*/
static int calculate_gas_optics(GasOptics_t * const lbl, /*Gas optics object.*/
                                Atmosphere_t const atm, /*Atmospheric state.*/
                                int const column, /*Column index.*/
                                int const offset, /*Additional index offset.*/
                                Optics_t * const gas, /*Optics due to molecular lines.*/
                                Optics_t * const rayleigh /*Optics due to rayleigh scattering.*/
                               )
{
    int i = column + offset;
    fp_t *level_pressure = &(atm.level_pressure[i*atm.num_levels]);
    fp_t *level_temperature = &(atm.level_temperature[i*atm.num_levels]);
    int j;
    for (j=0; j<atm.num_molecules; ++j)
    {
        fp_t const *ppmv = atm.ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_molecule_ppmv(lbl, atm.molecules[j], ppmv));
    }
    for (j=0; j<atm.num_cfcs; ++j)
    {
        fp_t const *ppmv = atm.cfc_ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_cfc_ppmv(lbl, atm.cfc[j].id, ppmv));
    }
    for (j=0; j<atm.num_cia_species; ++j)
    {
        fp_t const *ppmv = atm.cia_ppmv[j];
        ppmv = &(ppmv[i*atm.num_levels]);
        catch(set_cia_ppmv(lbl, atm.cia_species[j], ppmv));
    }
    catch(calculate_optical_depth(lbl, level_pressure, level_temperature, gas));
    catch(rayleigh_scattering(rayleigh, level_pressure));
    return GRTCODE_SUCCESS;
}


static int driver(Atmosphere_t const atm, /*Atmospheric state.*/
                  char const * const hitran_path, /*Path to HITRAN database ascii file.*/
                  char const * const solar_flux_path, /*Path to solar flux ascii file.*/
                  SpectralGrid_t const lw_grid, /*Longwave spectral grid.*/
                  SpectralGrid_t const sw_grid, /*Shortwave spectral grid.*/
                  Device_t const device[16], /*Device object.*/
                  Output_t * const output, /*Output object.*/
                  char const * const beta_path, /*Path to beta distribution input file.*/
                  char const * const ice_path, /*Path to ice cloud parameterization input file.*/
                  char const * const liquid_path, /*Path to liquid cloud parameterization input file.*/
                  int col_s, /*Column starting index.*/
                  int col_e, /*Column ending index.*/
                  int num_gpus)
{
    /*Create a grid that spans the entire wavenumber range.*/
    SpectralGrid_t grid;
    double w0 = lw_grid.w0 < sw_grid.w0 ? lw_grid.w0 : sw_grid.w0;
    double wn = sw_grid.wn > lw_grid.wn ? sw_grid.wn : lw_grid.wn;
    double res = sw_grid.dw > lw_grid.dw ? sw_grid.dw : lw_grid.dw;
/*
    catch(create_spectral_grid(&grid, w0, wn, res));
*/
    catch(create_spectral_grid(&grid, lw_grid.w0, lw_grid.wn, lw_grid.dw));

    /*Intialize gas optics objects.*/
    GasOptics_t lbl_lw[16];
/*
    int method = line_sweep;
*/
    int method = line_sample;
    int i;
    for (i=0; i<num_gpus; ++i)
    {
        catch(create_gas_optics(&(lbl_lw[i]), atm.num_levels, &lw_grid, &(device[i]),
                                hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
        catch(add_molecules(&(lbl_lw[i]), atm));
    }
/*
    GasOptics_t lbl_sw;
    catch(create_gas_optics(&lbl_sw, atm.num_levels, &sw_grid, &device,
                            hitran_path, atm.h2o_ctm, atm.o3_ctm, NULL, &method));
    catch(add_molecules(&lbl_sw, atm));
*/

    /*Initialize optics objects.*/
    Optics_t optics_lw_gas[16];
    for (i=0; i<num_gpus; ++i)
    {
        catch(create_optics(&(optics_lw_gas[i]), atm.num_layers, &lw_grid,
                            &(device[i])));
    }
/*
    Optics_t optics_sw_gas;
    catch(create_optics(&optics_sw_gas, atm.num_layers, &sw_grid, &device));
*/
    Optics_t optics_lw_rayleigh[16];
    for (i=0; i<num_gpus; ++i)
    {
        catch(create_optics(&(optics_lw_rayleigh[i]), atm.num_layers,
                            &lw_grid, &(device[i])));
    }
/*
    Optics_t optics_sw_rayleigh;
    catch(create_optics(&optics_sw_rayleigh, atm.num_layers, &sw_grid, &device));
*/
    Optics_t optics_lw_aerosol;
/*
    Optics_t optics_sw_aerosol;
*/
    if (!atm.clean)
    {
/*
        catch(create_optics(&optics_lw_aerosol, atm.num_layers, &lw_grid, &device));
        catch(create_optics(&optics_sw_aerosol, atm.num_layers, &sw_grid, &device));
*/
    }


    CloudOptics_t optics_cloud[16];
    if (!atm.clear)
    {
        /*Initialize clouds library.*/
        if (beta_path[0] == '\0' || ice_path[0] == '\0' || liquid_path[0] == '\0')
        {
            fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
            fprintf(stderr, "-beta-path, -ice-path, and -liquid-path args required when"
                            " running with clouds.");
            return EXIT_FAILURE;
        }
        for (i=0; i<num_gpus; ++i)
        {
            catch(create_cloud_optics(&(optics_cloud[i]), beta_path, liquid_path,
                                      ice_path, atm.num_layers, lw_grid, device[i]));
        }
    }


    /*Initialize total optics.*/
    Optics_t optics[16];
    for (i=0; i<num_gpus; ++i)
    {
        catch(create_optics(&(optics[i]), atm.num_layers, &grid, &(device[i])));
    }

    /*Initialize the solar fluxe object.*/
/*
    SolarFlux_t solar_flux;
    catch(create_solar_flux(&solar_flux, &grid, solar_flux_path));
*/

    /*Initialize solver objects.*/
    Longwave_t longwave[16];
    for (i=0; i<num_gpus; ++i)
    {
        catch(create_longwave(&(longwave[i]), atm.num_levels, &lw_grid, &(device[i])));
    }
/*
    Shortwave_t shortwave;
    catch(create_shortwave(&shortwave, atm.num_levels, &grid, &device));
*/

    /*Initialize other buffers.*/
/*  fp_t *albedo_dir = (fp_t *)malloc(sizeof(*albedo_dir)*grid.n);*/
    fp_t * surface_emissivity = (fp_t *)malloc(sizeof(*surface_emissivity)*lw_grid.n);
/*
    fp_t * lw_flux_up = (fp_t *)malloc(sizeof(*lw_flux_up)*atm.num_levels*lw_grid.n);
    fp_t * lw_flux_down = (fp_t *)malloc(sizeof(*lw_flux_down)*atm.num_levels*lw_grid.n);
*/
/*
    fp_t * sw_flux_up = (fp_t *)malloc(sizeof(*sw_flux_up)*atm.num_levels*grid.n);
    fp_t * sw_flux_down = (fp_t *)malloc(sizeof(*sw_flux_down)*atm.num_levels*grid.n);
*/

    /*Loop through the columns.*/
    int t;
    for (t=0; t<atm.num_times; ++t)
    {
        int offset = t*atm.num_columns;
        int i;

        /*Write out zeros for the fluxes to make combining the files easier.*/
        fp_t * buffer = (fp_t *)malloc(sizeof(*buffer)*atm.num_levels*lw_grid.n);
        for (i=0; i<atm.num_levels*lw_grid.n; ++i)
        {
            buffer[i] = 0.;
        }
        for (i=0; i<atm.num_columns; ++i)
        {
            write_output(output, RLUTCSAF, buffer, t, i);
            write_output(output, RLUSCSAF, buffer, t, i);
            write_output(output, RLDSCSAF, buffer, t, i);
            write_output(output, RLUTAF, buffer, t, i);
            write_output(output, RLUSAF, buffer, t, i);
            write_output(output, RLDSAF, buffer, t, i);
        }
        free(buffer);

        /*Run calculate on for columns this rank is responsible for.*/
#pragma omp parallel for if(device[0] != HOST_ONLY) num_threads(num_gpus) default(shared)
        for (i=col_s; i<col_e; ++i)
        {
            int id = 0;
#ifdef _OPENMP
            if (device[0] != HOST_ONLY)
            {
                id = omp_get_thread_num();
            }
#endif
            /*Thread specific buffers for output.*/
            fp_t * lw_flux_up = (fp_t *)malloc(sizeof(*lw_flux_up)*atm.num_levels*lw_grid.n);
            fp_t * lw_flux_down = (fp_t *)malloc(sizeof(*lw_flux_down)*atm.num_levels*lw_grid.n);
            fp_t * lw_flux_up_sum = NULL;
            fp_t * lw_flux_down_sum = NULL;
            if (!atm.clear)
            {
                lw_flux_up_sum = (fp_t *)malloc(sizeof(*lw_flux_up_sum)*atm.num_levels*lw_grid.n);
                lw_flux_down_sum = (fp_t *)malloc(sizeof(*lw_flux_down_sum)*atm.num_levels*lw_grid.n);
                int k;
                for (k=0; k<atm.num_levels*lw_grid.n; ++k)
                {
                    lw_flux_up_sum[k] = 0.;
                    lw_flux_down_sum[k] = 0.;
                }
            }

            /*Longwave gas optics.*/
            Optics_t optics_lw_total;
            calculate_gas_optics(&(lbl_lw[id]), atm, i, offset, &(optics_lw_gas[id]),
                                 &(optics_lw_rayleigh[id]));
            Optics_t * optics_array[4] = {&(optics_lw_gas[id]), &(optics_lw_rayleigh[id])};
            add_optics(optics_array, 2, &optics_lw_total);

            /*Longwave clear-clean-sky fluxes.*/
            fp_t surface_temperature = atm.surface_temperature[offset+i];
            fp_t * layer_temperature = &(atm.layer_temperature[(offset+i)*atm.num_layers]);
            fp_t * level_temperature = &(atm.level_temperature[(offset+i)*atm.num_levels]);
            int n = (offset + i)*atm.emissivity_grid_size;
            fp_t const * emissivity = &(atm.surface_emissivity[n]);
            interpolate_to_grid(lw_grid, atm.emissivity_grid, emissivity,
                                atm.emissivity_grid_size, surface_emissivity,
                                linear_sample, constant_extrapolation);
            calculate_lw_fluxes(&(longwave[id]), &optics_lw_total, surface_temperature,
                                layer_temperature, level_temperature,
                                surface_emissivity, lw_flux_up, lw_flux_down);
/*
            write_output(output, TAU, optics_lw_total.tau, t, i);
            fp_t tauz[lw_grid.n];
            for (n=0; n<lw_grid.n; ++n)
            {
                tauz[n] = 0.;
            }
            for (n=0; n<atm.num_layers; ++n)
            {
                int j;
                for (j=0; j<lw_grid.n; ++j)
                {
                    tauz[j] += optics_lw_total.tau[n*lw_grid.n + j];
                }
            }
            write_output(output, TAUZ, tauz, t, i);
*/
            destroy_optics(&optics_lw_total);
/*
            write_output(output, RLUCSAF, lw_flux_up, t, i);
            write_output(output, RLDCSAF, lw_flux_down, t, i);
*/
            n = lw_grid.n*(atm.num_levels - 1);
#pragma omp critical
            {
                write_output(output, RLUTCSAF, lw_flux_up, t, i);
                write_output(output, RLUSCSAF, &(lw_flux_up[n]), t, i);
                write_output(output, RLDSCSAF, &(lw_flux_up[n]), t, i);
            }

            if (!atm.clear)
            {
                /*Calculate overlap for the column.*/
                fp_t const * pressure = &(atm.layer_pressure[(offset+i)*atm.num_layers]);
                fp_t const * cloud_fraction = &(atm.cloud_fraction[(offset+i)*atm.num_layers]);
                fp_t const * liquid_content = &(atm.liquid_water_content[(offset+i)*atm.num_layers]);
                fp_t const * ice_content = &(atm.ice_water_content[(offset+i)*atm.num_layers]);
                fp_t const * thickness = &(atm.layer_thickness[(offset+i)*atm.num_layers]);

                /*Loop over the subcolumns.*/
                int j;
                int const num_subcolumns = 10;
                for (j=0; j<num_subcolumns; ++j)
                {
                    calculate_cloud_optics(&(optics_cloud[id]), atm.num_layers, pressure,
                                           layer_temperature, thickness, cloud_fraction,
                                           liquid_content, ice_content);

                    /*Add to gas optics.*/
                    optics_array[2] = &(optics_cloud[id].liquid_optics);
                    optics_array[3] = &(optics_cloud[id].ice_optics);
                    add_optics(optics_array, 4, &optics_lw_total);

                    /*Longwave aerosol-free fluxes.*/
                    calculate_lw_fluxes(&(longwave[id]), &optics_lw_total, surface_temperature,
                                        layer_temperature, level_temperature,
                                        surface_emissivity, lw_flux_up, lw_flux_down);
                    destroy_optics(&optics_lw_total);

                    /*Add the fluxes to the running subcolumn sums.*/
                    int k;
                    for (k=0; k<atm.num_levels*lw_grid.n; ++k)
                    {
                        lw_flux_up_sum[k] += lw_flux_up[k];
                        lw_flux_down_sum[k] += lw_flux_down[k];
                    }
                }
                /*Divide running subcolumn flux sums by the number of subcolumns.*/
                int k;
                for (k=0; k<atm.num_levels*lw_grid.n; ++k)
                {
                    lw_flux_up[k] = lw_flux_up_sum[k]/((fp_t)num_subcolumns);
                    lw_flux_down[k] = lw_flux_down_sum[k]/((fp_t)num_subcolumns);
                }
                /*Write out the fluxes.*/
#pragma omp critical
                {
                    write_output(output, RLUTAF, lw_flux_up, t, i);
                    write_output(output, RLUSAF, &lw_flux_up[n], t, i);
                    write_output(output, RLDSAF, &lw_flux_up[n], t, i);
                }
            }

            free(lw_flux_up);
            free(lw_flux_down);
            if (!atm.clear)
            {
                free(lw_flux_up_sum);
                free(lw_flux_down_sum);
            }
            continue;

#ifdef FOOBAR
            if (!atm.clean)
            {
                /*Longwave aerosol optics.*/
                catch(calculate_aerosol_optics(atm, i, offset, &optics_lw_aerosol));
                optics_array[2] = &optics_lw_aerosol;
                catch(add_optics(optics_array, 3, &optics_lw_total));

                /*Longwave clear-sky fluxes.*/
/*
                catch(calculate_lw_fluxes(&longwave, &optics_lw_total, surface_temperature,
                                          layer_temperature, level_temperature,
                                          surface_emissivity, lw_flux_up, lw_flux_down));
                catch(destroy_optics(&optics_lw_total));
                write_output(output, RLUCS, lw_flux_up, t, i);
                write_output(output, RLDCS, lw_flux_down, t, i);
*/
            }

            write_output(output, PLEV, &(atm.level_pressure[(offset + i)*atm.num_levels]), t, i);
            write_output(output, TLEV, &(atm.level_temperature[(offset + i)*atm.num_levels]), t, i);
            write_output(output, TLAY, &(atm.layer_temperature[(offset + i)*atm.num_layers]), t, i);
            write_output(output, TS, &(atm.surface_temperature[offset + i]), t, i);
            int j;
            for (j=0; j<atm.num_molecules; ++j)
            {
                fp_t *ppmv = atm.ppmv[j];
                write_output(output, (VarId_t)((int)H2OVMR + j),
                              &(ppmv[(offset + i)*atm.num_levels]), t, i);
            }
            continue;
#endif

            /*Shortwave.*/
#ifdef SHORTWAVE
            fp_t const zen_dir = atm.solar_zenith_angle[offset + i];
            if (zen_dir > 0.)
            {
                Optics_t optics_sw_total;
                catch(calculate_optics(&lbl_sw, atm, i, offset, &optics_sw_gas,
                                       &optics_sw_rayleigh, &optics_sw_aerosol,
                                       &optics_sw_cloud, &optics_sw_total));
                /*Combine the longwave and shortwave optics.*/
                catch(sample_optics(&optics, &optics_lw_total, &lw_grid.w0, &lw_grid.wn));
                catch(sample_optics(&optics, &optics_sw_total, &sw_grid.w0, &sw_grid.wn));
                fp_t const zen_dif = 0.5;
                int n = (offset + i)*atm.albedo_grid_size;
                catch(interpolate_to_grid(grid, atm.albedo_grid, &(atm.surface_albedo[n]),
                                          atm.albedo_grid_size, albedo_dir, linear_sample,
                                          constant_extrapolation));
                fp_t *albedo_dif = albedo_dir;
                catch(calculate_sw_fluxes(&shortwave, &optics, zen_dir, zen_dif,
                                          albedo_dir, albedo_dif, atm.total_solar_irradiance[offset+i],
                                          solar_flux.incident_flux, sw_flux_up, sw_flux_down));
                catch(destroy_optics(&optics_sw_total));

                /*Integrate fluxes and write them to the output file.*/
                for (j=0; j<atm.num_levels; ++j)
                {
                    integrate(&(sw_flux_up[j*grid.n]), grid.n, grid.dw, &(flux_up_total[j]));
                    integrate(&(sw_flux_down[j*grid.n]), grid.n, grid.dw, &(flux_down_total[j]));
                }
/*
                write_output(output, RSU, flux_up_total, t, i);
                write_output(output, RSD, flux_down_total, t, i);
*/
            }
#endif
        }
    }

    /*Release memory.*/
/*  free(albedo_dir);*/
    free(surface_emissivity);
/*
    free(lw_flux_up);
    free(lw_flux_down);
*/
/*
    free(sw_flux_up);
    free(sw_flux_down);
    catch(destroy_shortwave(&shortwave));
*/
    for (i=0; i<num_gpus; ++i)
    {
        catch(destroy_longwave(&(longwave[i])));
    }
/*
    catch(destroy_solar_flux(&solar_flux));
*/
    for (i=0; i<num_gpus; ++i)
    {
        catch(destroy_optics(&(optics_lw_gas[i])));
    }
/*
    catch(destroy_optics(&optics_sw_gas));
*/
    for (i=0; i<num_gpus; ++i)
    {
        catch(destroy_optics(&(optics_lw_rayleigh[i])));
    }
/*
    catch(destroy_optics(&optics_sw_rayleigh));
*/
    if (!atm.clean)
    {
/*
        catch(destroy_optics(&optics_lw_aerosol));
        catch(destroy_optics(&optics_sw_aerosol));
*/
    }
    if (!atm.clear)
    {
        for (i=0; i<num_gpus; ++i)
        {
            catch(destroy_cloud_optics(&(optics_cloud[i])));
        }
    }
    for (i=0; i<num_gpus; ++i)
    {
        catch(destroy_optics(&(optics[i])));
        catch(destroy_gas_optics(&(lbl_lw[i])));
    }
/*
    catch(destroy_gas_optics(&lbl_sw));
*/
    return GRTCODE_SUCCESS;
}


/*Main driver program.  When linking an executable, the user must provide an object file that
  includes implementations for the following defined in driver.h:
   - struct output;
   - void close_flux_file(Output_t *output);
   - Atmosphere_t create_atmosphere(Parser_t parser);
   - void destroy_atmosphere(Atmosphere_t *atm);
   - Output_t create_flux_file(char *path, Atmosphere_t *atm);
   - void write_output(Output_t *output, enum varid id, fp_t *data, int time, int column);*/
int main(int argc, char **argv)
{
    /*Add/parse command line arguments.*/
    char *description = "Calculates radiative fluxes using the line-by-line method.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "hitran_file", NULL, "HITRAN database file.", NULL);
    add_argument(&parser, "solar_flux", NULL, "Solar flux file.", NULL);
    int one = 1;
    add_argument(&parser, "-beta-path", NULL, "Path to beta distribution input file.", &one);
    add_argument(&parser, "-c", "--line-cutoff", "Cutoff [1/cm] from line center.", &one);
    add_argument(&parser, "-ice-path", NULL, "Path to ice cloud parameterization input file.", &one);
    add_argument(&parser, "-liquid-path", NULL, "Path to liquid cloud parameterization input file.", &one);
    add_argument(&parser, "-o", NULL, "Name of output file.", &one);
    add_argument(&parser, "-r-lw", "--lw-resolution", "Longwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-r-sw", "--sw-resolution", "Shortwave spectral resolution [1/cm].", &one);
    add_argument(&parser, "-s", "--solver", "Shortwave solver.", &one);
    add_argument(&parser, "-v", "--verbose", "Increase verbosity.", NULL);
    add_argument(&parser, "-w-lw", "--lw-lower-bound", "Longwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-w-sw", "--sw-lower-bound", "Shortwave spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-W-lw", "--lw-upper-bound", "Longwave spectral upper bound [1/cm].", &one);
    add_argument(&parser, "-W-sw", "--sw-upper-bound", "Shortwave spectral upper bound [1/cm].", &one);

    /*Get the atmospheric data.*/
    Atmosphere_t atm = create_atmosphere(&parser);

    /*Start up MPI.*/
    int rank;
    int col_s;
    int col_e;
    distribute_init(atm.num_columns, &rank, &col_s, &col_e);

    /*Set GRTCODE verbosity.*/
    int verbosity_level = get_argument(parser, "-v", NULL) ? GRTCODE_INFO : GRTCODE_WARN;
    grtcode_set_verbosity(verbosity_level);

    /*Get paths to required input files.*/
    char hitran_path[valuelen];
    get_argument(parser, "hitran_file", hitran_path);
    char solar_flux_path[valuelen];
    get_argument(parser, "solar_flux", solar_flux_path);

    /*Create a spectral grids.*/
    char buffer[valuelen];
    double w0 = get_argument(parser, "-w-lw", buffer) ? atof(buffer) : 1.;
    double wn = get_argument(parser, "-W-lw", buffer) ? atof(buffer) : 3250.;
    double dw = get_argument(parser, "-r-lw", buffer) ? atof(buffer) : 0.1;
    SpectralGrid_t lw_grid;
    catch(create_spectral_grid(&lw_grid, w0, wn, dw));
    dw = get_argument(parser, "-r-sw", buffer) ? atof(buffer) : 1.;
    w0 = get_argument(parser, "-w-sw", buffer) ? atof(buffer) : dw + wn;
    wn = get_argument(parser, "-W-sw", buffer) ? atof(buffer) : 50000.;
    SpectralGrid_t sw_grid;
    catch(create_spectral_grid(&sw_grid, w0, wn, dw));

    /*Set the device to run on.*/
    int num_gpus;
    catch(get_num_gpus(&num_gpus, 0));
    Device_t device[16];
    if (num_gpus == 0)
    {
        /*Run on the host-only.*/
        num_gpus = 1;
        catch(create_device(&(device[0]), NULL));
    }
    else
    {
        int i;
        for (i=0; i<num_gpus; ++i)
        {
            catch(create_device(&(device[i]), &i));
        }
        /*Set one openmp thread per GPU.*/
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(num_gpus);
#endif
    }

    /*Initialize the output file.*/
    char path[valuelen];
    if (get_argument(parser, "-o", buffer))
    {
        snprintf(path, valuelen, "%d.%s", rank, buffer);
    }
    else
    {
        snprintf(path, valuelen, "%d.%s", rank, "output.nc");
    }
    Output_t *output;
    create_flux_file(&output, path, &atm, &lw_grid, &sw_grid);

    /*Get cloud parameterization inputs.*/
    char beta_path[valuelen];
    if (!get_argument(parser, "-beta-path", beta_path))
    {
        beta_path[0] = '\0';
    }
    char ice_path[valuelen];
    if (!get_argument(parser, "-ice-path", ice_path))
    {
        ice_path[0] = '\0';
    }
    char liquid_path[valuelen];
    if (!get_argument(parser, "-liquid-path", liquid_path))
    {
        liquid_path[0] = '\0';
    }

    /*Calculate the fluxes over all the columns.*/
    catch(driver(atm, hitran_path, solar_flux_path, lw_grid, sw_grid, device, output,
                 beta_path, ice_path, liquid_path, col_s, col_e, num_gpus));

    /*Clean up.*/
    close_flux_file(output);
    destroy_atmosphere(&atm);
    destroy_parser(&parser);

    /*Finalize MPI.*/
    distribute_final();

    return EXIT_SUCCESS;
}
