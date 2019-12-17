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
#include "atmosphere.h"
#include "gas_optics.h"
#include "grtcode_utilities.h"
#include "longwave.h"
#include "rayleigh.h"
#include "shortwave.h"
#include "solar_flux.h"


#define catch(e) { \
    int e_ = e; \
    if (e_ != GRTCODE_SUCCESS) { \
        fprintf(stderr, "[%s, %d] Error:\n", __FILE__, __LINE__); \
        char b_[1024]; \
        grtcode_errstr(e_, b_, 1024); \
        fprintf(stderr, "%s", b_); \
        return EXIT_FAILURE; \
    }}
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define MAX_NUM_MOLECULES 7
#define MAX_NUM_CFCS 21
#define MAX_NUM_CIAS 3


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


/*Set a molecular species as active.*/
static void activate_molecule(Parser_t const parser, /*Parser object.*/
                              char const * const arg, /*Arg to check for.*/
                              int * const array, /*Array.*/
                              int const tag, /*Value to store in array.*/
                              int * spot, /*Current spot in array.*/
                              int const max_size /*Size of input array.*/
                             )
{
    char buffer[valuelen];
    if (get_argument(parser, arg, buffer))
    {
        if (*spot >= max_size)
        {
            fprintf(stderr, "Array is too small, increase size.\n");
            exit(EXIT_FAILURE);
        }
        array[*spot] = tag;
        *spot += 1;
    }
    return;
}


/*Set a CFC as active.*/
static void activate_cfc(Parser_t const parser, /*Parser object.*/
                         char const * const arg, /*Arg to check for.*/
                         char const * const arg_eq, /*Equivalent arg.*/
                         Cfc_t * const array, /*Array of Cfc_t objects.*/
                         int const id, /*CFC id to store.*/
                         int * spot, /*Current spot in array.*/
                         int const max_size)
{
    int found;
    int use_eq = 0;
    char buffer[valuelen];
    found = get_argument(parser, arg, buffer);
    if (!found && arg_eq != NULL)
    {
        use_eq = get_argument(parser, arg_eq, buffer);
    }
    if (found || use_eq)
    {
        if (*spot >= max_size)
        {
            fprintf(stderr, "Array is too small, increase size.\n");
            exit(EXIT_FAILURE);
        }
        array[*spot].id = id;
        snprintf(array[*spot].path, valuelen, "%s", buffer);
        array[*spot].use_equivalent_ppmv = use_eq;
        *spot += 1;
    }
    return;
}


/*Set collision-induced absorption as active.*/
static void activate_cia(Parser_t const parser, /*Parser object.*/
                         char const * const arg, /*Argument to check for.*/
                         int * const array, /*Array of species ids.*/
                         int const id1, /*Id of species.*/
                         int const id2, /*Id of species.*/
                         int * spot, /*Current spot in input array of species ids.*/
                         int const max_size, /*Size of input array of species ids.*/
                         char **values, /*Array to store argument values.*/
                         int * values_spot, /*Current spot in input are of values.*/
                         int * combos /*Store id combinations.*/
                        )
{
    char buffer[valuelen];
    if (get_argument(parser, arg, buffer))
    {
        if (*values_spot >= max_size)
        {
            fprintf(stderr, "Input values array is too small, increase size.\n");
            exit(EXIT_FAILURE);
        }
        snprintf(values[*values_spot], valuelen, "%s", buffer);
        int index = 2*(*values_spot);
        combos[index] = id1;
        combos[index+1] = id2;
        *values_spot += 1;
        int ids[2] = {id1, id2};
        int i;
        for (i=0; i<2; ++i)
        {
            if (*spot >= max_size)
            {
                fprintf(stderr, "Input array is too small, increase size.\n");
                exit(EXIT_FAILURE);
            }
            int found = 0;
            int j;
            for (j=0; j<*spot; ++j)
            {
                if (ids[i] == array[j])
                {
                    found = 1;
                    break;
                }
            }
            if (!found)
            {
                array[*spot] = ids[i];
                *spot += 1;
            }
        }
    }
    return;
}


/*Calculate the radiative fluxes for the CMIP6 RFMIP-IRF test cases.*/
int main(int argc, char **argv)
{
    /*Add/parse command line arguments.*/
    char *description = "Calculates the radiative fluxes for the CMIP6 RFMIP-IRF test"
                        " cases.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "input_file", NULL, "Input data file.", NULL);
    add_argument(&parser, "experiment", NULL, "Experiment number.", NULL);
    add_argument(&parser, "hitran_file", NULL, "HITRAN database file.", NULL);
    add_argument(&parser, "solar_flux", NULL, "Solar flux CSV file.", NULL);
    int one = 1;
    add_argument(&parser, "-CCl4", NULL, "CSV file with CCl4 cross sections.", &one);
    add_argument(&parser, "-C2F6", NULL, "CSV file with C2F6 cross sections.", &one);
    add_argument(&parser, "-CF4", NULL, "CSV file with CF4 cross sections.", &one);
    add_argument(&parser, "-CFC-11", NULL, "CSV file with CFC-11 cross sections.", &one);
    add_argument(&parser, "-CFC-11-eq", NULL, "CSV file with CFC-11 cross sections.", &one);
    add_argument(&parser, "-CFC-12", NULL, "CSV file with CFC-12 cross sections.", &one);
    add_argument(&parser, "-CFC-12-eq", NULL, "CSV file with CFC-12 cross sections.", &one);
    add_argument(&parser, "-CFC-113", NULL, "CSV file with CFC-113 cross sections.", &one);
    add_argument(&parser, "-CFC-114", NULL, "CSV file with CFC-114 cross sections.", &one);
    add_argument(&parser, "-CFC-115", NULL, "CSV file with CFC-115 cross sections.", &one);
    add_argument(&parser, "-CH2Cl2", NULL, "CSV file with CH2Cl2 cross sections.", &one);
    add_argument(&parser, "-CH4", NULL, "Include CH4.", NULL);
    add_argument(&parser, "-CO", NULL, "Include CO.", NULL);
    add_argument(&parser, "-CO2", NULL, "Include CO2.", NULL);
    add_argument(&parser, "-H2O", NULL, "Include H2O.", NULL);
    add_argument(&parser, "-HCFC-22", NULL, "CSV file with HCFC-22 cross sections.", &one);
    add_argument(&parser, "-HCFC-141b", NULL, "CSV file with HCFC-141b cross sections.", &one);
    add_argument(&parser, "-HCFC-142b", NULL, "CSV file with HCFC-142b cross sections.", &one);
    add_argument(&parser, "-HFC-23", NULL, "CSV file with HFC-23 cross sections.", &one);
    add_argument(&parser, "-HFC-125", NULL, "CSV file with HFC-125 cross sections.", &one);
    add_argument(&parser, "-HFC-134a", NULL, "CSV file with HFC-134a cross sections.", &one);
    add_argument(&parser, "-HFC-134a-eq", NULL, "CSV file with HFC-134a cross sections.", &one);
    add_argument(&parser, "-HFC-143a", NULL, "CSV file with HFC-143a cross sections.", &one);
    add_argument(&parser, "-HFC-152a", NULL, "CSV file with HFC-152a cross sections.", &one);
    add_argument(&parser, "-HFC-227ea", NULL, "CSV file with HFC-227ea cross sections.", &one);
    add_argument(&parser, "-HFC-245fa", NULL, "CSV file with HFC-245fa cross sections.", &one);
    add_argument(&parser, "-N2-N2", NULL, "CSV file with N2-N2 collison cross sections", &one);
    add_argument(&parser, "-N2O", NULL, "Include N2O.", NULL);
    add_argument(&parser, "-NF3", NULL, "CSV file with NF3 cross sections.", &one);
    add_argument(&parser, "-O2", NULL, "Include O2.", NULL);
    add_argument(&parser, "-O2-N2", NULL, "CSV file with O2-N2 collison cross sections", &one);
    add_argument(&parser, "-O2-O2", NULL, "CSV file with O2-O2 collison cross sections", &one);
    add_argument(&parser, "-O3", NULL, "Include O3.", NULL);
    add_argument(&parser, "-SF6", NULL, "CSV file with SF6 cross sections.", &one);
    add_argument(&parser, "-c", "--line-cutoff", "Cutoff [1/cm] from line center.", &one);
    add_argument(&parser, "-d", "--device", "GPU id", &one);
    add_argument(&parser, "-h2o-ctm", NULL, "Directory containing H2O continuum files", &one);
    add_argument(&parser, "-o", NULL, "Name of output file.", &one);
    add_argument(&parser, "-o3-ctm", NULL, "Directory containing O3 continuum files", &one);
    add_argument(&parser, "-r", "--spectral-resolution", "Spectral resolution [1/cm].", &one);
    add_argument(&parser, "-v", "--verbose", "Increase verbosity.", NULL);
    add_argument(&parser, "-w", "--spectral-lower-bound", "Spectral lower bound [1/cm].", &one);
    add_argument(&parser, "-W", "--spectral-upper-bound", "Spectral upper bound [1/cm].", &one);
    add_argument(&parser, "-x", "--column-lower-bound", "Starting column index.", &one);
    add_argument(&parser, "-X", "--column-upper-bound", "Ending column index.", &one);
    add_argument(&parser, "-z", "--level-lower-bound", "Starting level index.", &one);
    add_argument(&parser, "-Z", "--level-upper-bound", "Ending level index.", &one);
    parse_args(parser);

    /*Set verbosity.*/
    char buffer[valuelen];
    if (get_argument(parser, "-v", NULL))
    {
        grtcode_set_verbosity(GRTCODE_INFO);
    }
    else
    {
        grtcode_set_verbosity(GRTCODE_WARN);
    }

    /*Set device.*/
    Device_t device;
    if (get_argument(parser, "-d", buffer))
    {
        int d = atoi(buffer);
        catch(create_device(&device, &d));
    }
    else
    {
        catch(create_device(&device, NULL));
    }

    /*Create a spectral grid.*/
    double w0 = 1.;
    if (get_argument(parser, "-w", buffer))
    {
        w0 = atof(buffer);
    }
    double wn = 50000.;
    if (get_argument(parser, "-W", buffer))
    {
        wn = atof(buffer);
    }
    double dw = 0.1;
    if (get_argument(parser, "-r", buffer))
    {
        dw = atof(buffer);
    }
    SpectralGrid_t grid;
    catch(create_spectral_grid(&grid, w0, wn, dw));

    /*Determine which molecules to use.*/
    int molecules[MAX_NUM_MOLECULES];
    int num_molecules = 0;
    activate_molecule(parser, "-CH4", molecules, CH4, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-CO", molecules, CO, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-CO2", molecules, CO2, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-H2O", molecules, H2O, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-N2O", molecules, N2O, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-O2", molecules, O2, &num_molecules, MAX_NUM_MOLECULES);
    activate_molecule(parser, "-O3", molecules, O3, &num_molecules, MAX_NUM_MOLECULES);

    /*Determine which CFCs to use.  "Equivalent" concentrations will override
      non-equivalent ones if both are flags are used.*/
    Cfc_t cfc[MAX_NUM_CFCS];
    int i;
    for (i=0; i<MAX_NUM_CFCS; ++i)
    {
        cfc[i].path = (char *)malloc(sizeof(*(cfc[i].path))*valuelen);
    }
    int num_cfcs = 0;
    activate_cfc(parser, "-CCl4", NULL, cfc, CCl4, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-C2F6", NULL, cfc, C2F6, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CF4", NULL, cfc, CF4, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CH2Cl2", NULL, cfc, CH2Cl2, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CFC-11", "-CFC-11-eq", cfc, CFC11, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CFC-12", "-CFC-12-eq", cfc, CFC12, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CFC-113", NULL, cfc, CFC113, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CFC-114", NULL, cfc, CFC114, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-CFC-115", NULL, cfc, CFC115, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HCFC-22", NULL, cfc, HCFC22, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HCFC-141b", NULL, cfc, HCFC141b, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HCFC-142b", NULL, cfc, HCFC142b, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-23", NULL, cfc, HFC23, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-125", NULL, cfc, HFC125, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-134a", "-HFC-134a-eq", cfc, HFC134a, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-143a", NULL, cfc, HFC143a, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-152a", NULL, cfc, HFC152a, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-227ea", NULL, cfc, HFC227ea, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-HFC-245fa", NULL, cfc, HFC245fa, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-NF3", NULL, cfc, NF3, &num_cfcs, MAX_NUM_CFCS);
    activate_cfc(parser, "-SF6", NULL, cfc, SF6, &num_cfcs, MAX_NUM_CFCS);

    /*Determine which collision-induced absorption spectra to include.*/
    int cia_species[MAX_NUM_CIAS];
    char *cia_path[MAX_NUM_CIAS];
    int cia_combos[2*MAX_NUM_CIAS];
    for (i=0; i<MAX_NUM_CIAS; ++i)
    {
        cia_path[i] = (char *)malloc(sizeof(*(cia_path[i]))*valuelen);
    }
    int num_cia_species = 0;
    int num_cias = 0;
    activate_cia(parser, "-N2-N2", cia_species, CIA_N2, CIA_N2, &num_cia_species,
                 MAX_NUM_CIAS, cia_path, &num_cias, cia_combos);
    activate_cia(parser, "-O2-O2", cia_species, CIA_O2, CIA_O2, &num_cia_species,
                 MAX_NUM_CIAS, cia_path, &num_cias, cia_combos);
    activate_cia(parser, "-O2-N2", cia_species, CIA_O2, CIA_N2, &num_cia_species,
                 MAX_NUM_CIAS, cia_path, &num_cias, cia_combos);

    /*Read in the atmospheric input data.*/
    Atmosphere_t atm;
    atm.num_lw_wavenumber = grid.n;
    atm.num_sw_wavenumber = grid.n;
    atm.x = 0;
    if (get_argument(parser, "-x", buffer))
    {
        atm.x = atoi(buffer);
    }
    atm.num_columns = 100;
    if (get_argument(parser, "-X", buffer))
    {
        atm.num_columns = atoi(buffer) - atm.x + 1;
    }
    atm.z = 0;
    if (get_argument(parser, "-z", buffer))
    {
        atm.z = atoi(buffer);
    }
    atm.num_levels = 61;
    if (get_argument(parser, "-Z", buffer))
    {
        atm.num_levels = atoi(buffer) - atm.z + 1;
    }
    atm.num_layers = atm.num_levels - 1;
    get_argument(parser, "experiment", buffer);
    int experiment = atoi(buffer);
    get_argument(parser, "input_file", buffer);
    create_atmosphere(&atm, buffer, experiment, molecules, num_molecules, cfc,
                      num_cfcs, cia_species, num_cia_species);

    /*Read in the incident solar flux.*/
    SolarFlux_t solar_flux;
    get_argument(parser, "solar_flux", buffer);
    catch(create_solar_flux(&solar_flux, &grid, buffer));

    /*Initialize a molecular lines object.*/
    char hitran_path[valuelen];
    get_argument(parser, "hitran_file", hitran_path);
    char h2o_ctm[valuelen];
    if (!get_argument(parser, "-h2o-ctm", h2o_ctm))
    {
        snprintf(h2o_ctm, valuelen, "%s", "none");
    }
    char o3_ctm[valuelen];
    if (!get_argument(parser, "-o3-ctm", o3_ctm))
    {
        snprintf(o3_ctm, valuelen, "%s", "none");
    }
    GasOptics_t lbl;
    int method = line_sweep;
    catch(create_gas_optics(&lbl, atm.num_levels, &grid, &device,
                            hitran_path, h2o_ctm, o3_ctm, NULL, &method));

    /*Add molecules, CFCs, and collision-induced absoprtion.*/
    for (i=0; i<num_molecules; ++i)
    {
        catch(add_molecule(&lbl, molecules[i], NULL, NULL));
    }
    for (i=0; i<num_cfcs; ++i)
    {
        catch(add_cfc(&lbl, cfc[i].id, cfc[i].path));
    }
    for (i=0; i<num_cias; ++i)
    {
        int index = 2*i;
        catch(add_cia(&lbl, cia_combos[index], cia_combos[index+1],
                      cia_path[i]));
    }

    /*Initialize an optics object.*/
    Optics_t optics_ghgs;
    catch(create_optics(&optics_ghgs, atm.num_layers, &grid, &device));
    Optics_t optics_rayleigh;
    catch(create_optics(&optics_rayleigh, atm.num_layers, &grid, &device));

    /*Initialize a longwave object.*/
    Longwave_t longwave;
    catch(create_longwave(&longwave, atm.num_levels, &grid, &device));

    /*Initialize a shortwave object.*/
    Shortwave_t shortwave;
    catch(create_shortwave(&shortwave, atm.num_levels, &grid, &device));

    /*Initialize the output file.*/
    if (!get_argument(parser, "-o", buffer))
    {
        snprintf(buffer, valuelen, "%s", "rfmip-irf.output.nc");
    }
    Output_t output = create_flux_file(buffer, &atm);

    /*Loop through the columns.*/
    fp_t *flux_up = (fp_t *)malloc(sizeof(*flux_up)*atm.num_levels*grid.n);
    fp_t *flux_down = (fp_t *)malloc(sizeof(*flux_down)*atm.num_levels*grid.n);
    for (i=0; i<atm.num_columns; ++i)
    {
        /*Calculate molecular spectra.*/
        fp_t *level_pressure = &(atm.level_pressure[i*atm.num_levels]);
        fp_t *level_temperature = &(atm.level_temperature[i*atm.num_levels]);
        int j;
        for (j=0; j<num_molecules; ++j)
        {
            fp_t *ppmv = atm.ppmv[j];
            ppmv = &(ppmv[i*atm.num_levels]);
            catch(set_molecule_ppmv(&lbl, molecules[j], ppmv));
        }
        for (j=0; j<num_cfcs; ++j)
        {
            fp_t *ppmv = atm.cfc_ppmv[j];
            ppmv = &(ppmv[i*atm.num_levels]);
            catch(set_cfc_ppmv(&lbl, cfc[j].id, ppmv));
        }
        for (j=0; j<num_cia_species; ++j)
        {
            fp_t *ppmv = atm.cia_ppmv[j];
            ppmv = &(ppmv[i*atm.num_levels]);
            catch(set_cia_ppmv(&lbl, cia_species[j], ppmv));
        }
        catch(calculate_optical_depth(&lbl, level_pressure, level_temperature, &optics_ghgs));

        /*Calculate longwave fluxes.*/
        fp_t surface_temperature = atm.surface_temperature[i];
        fp_t *layer_temperature = &(atm.layer_temperature[i*atm.num_layers]);
        fp_t *surface_emissivity = &(atm.surface_emissivity[i*atm.num_lw_wavenumber]);
        double lw_solver_w0 = 1. < grid.w0 ? grid.w0 : 1.;
        double lw_solver_wn = 3250. > grid.wn ? grid.wn : 3250.;
        SpectralGrid_t lw_solver_grid;
        catch(create_spectral_grid(&lw_solver_grid, lw_solver_w0, lw_solver_wn, grid.dw));
        catch(calculate_lw_fluxes(&longwave, &optics_ghgs, surface_temperature,
                                  layer_temperature, level_temperature,
                                  surface_emissivity, flux_up, flux_down,
                                  &(lw_solver_grid.w0), &(lw_solver_grid.wn)));

        /*Integrate fluxes and write them to the output file.*/
        fp_t flux_up_total[atm.num_levels];
        fp_t flux_down_total[atm.num_levels];
        for (j=0; j<atm.num_levels; ++j)
        {
            integrate(&(flux_up[j*lw_solver_grid.n]), lw_solver_grid.n, lw_solver_grid.dw,
                      &(flux_up_total[j]));
            integrate(&(flux_down[j*lw_solver_grid.n]), lw_solver_grid.n, lw_solver_grid.dw,
                      &(flux_down_total[j]));
        }
        size_t start[3] = {0, i, 0};
        size_t count[3] = {1, 1, atm.num_levels};
        write_output(&output, RLU, flux_up_total, start, count);
        write_output(&output, RLD, flux_down_total, start, count);

        fp_t const zen_dir = atm.solar_zenith_angle[i];
        if (zen_dir > 0.)
        {
            /*Calculate the optical properities of a column.*/
            catch(rayleigh_scattering(&optics_rayleigh, level_pressure));

            /*Calculate the combined optical properties.*/
            Optics_t const * const optics_mech[2] = {&optics_ghgs, &optics_rayleigh};
            Optics_t optics_combined;
            catch(add_optics(optics_mech, 2, &optics_combined));

            /*Calculate shortwave fluxes.*/
            fp_t const zen_dif = 0.5;
            fp_t *albedo_dir = &(atm.surface_albedo[i*atm.num_sw_wavenumber]);
            fp_t *albedo_dif = albedo_dir;
            catch(calculate_sw_fluxes(&shortwave, &optics_combined, zen_dir, zen_dif,
                                      albedo_dir, albedo_dif, atm.total_solar_irradiance[i],
                                      solar_flux.incident_flux, flux_up, flux_down));
            catch(destroy_optics(&optics_combined));

            /*Integrate fluxes and write them to the output file.*/
            for (j=0; j<atm.num_levels; ++j)
            {
                integrate(&(flux_up[j*grid.n]), grid.n, grid.dw, &(flux_up_total[j]));
                integrate(&(flux_down[j*grid.n]), grid.n, grid.dw, &(flux_down_total[j]));
            }
            write_output(&output, RSU, flux_up_total, start, count);
            write_output(&output, RSD, flux_down_total, start, count);
        }
    }

    /*Clean up.*/
    close_flux_file(&output);
    free(flux_up);
    free(flux_down);
    catch(destroy_shortwave(&shortwave));
    catch(destroy_longwave(&longwave));
    catch(destroy_solar_flux(&solar_flux));
    catch(destroy_optics(&optics_ghgs));
    catch(destroy_optics(&optics_rayleigh));
    catch(destroy_gas_optics(&lbl));
    destroy_atmosphere(&atm);
    destroy_parser(&parser);
    for (i=0; i<MAX_NUM_CFCS; ++i)
    {
        free(cfc[i].path);
    }
    for (i=0; i<MAX_NUM_CIAS; ++i)
    {
        free(cia_path[i]);
    }
    return EXIT_SUCCESS;
}
