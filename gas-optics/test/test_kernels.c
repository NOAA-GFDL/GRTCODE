#include <stdint.h>
#include <string.h>
#include <curtis_godson.h>
#include "grtcode_utilities.h"
#include "kernels.h"
#include "kernel_utils.h"
#include "line_shape.h"
#include "molecules.h"
#include "RFM_voigt.h"
#include "spectral_bin.h"
#include "spectral_bin-internal.h"
#include "test_harness.h"
#include "tips2017.h"


#define NUM_TESTS 15
#define NUM_ISOTOPOLOGUES 9
#define NUM_LAYERS 5
#define NUM_LEVELS 6
#define NUM_TRANSITIONS 5
#define NUM_GRID_POINTS 50


/*Helper structure that mimics an atmospheric column*/
struct Atmosphere
{
    fp_t level_pressure[NUM_LEVELS];
    fp_t level_xh2o[NUM_LEVELS];
    fp_t pressure[NUM_LAYERS];
    fp_t temperature[NUM_LAYERS];
    fp_t xh2o[NUM_LAYERS];
    fp_t number_density[NUM_LAYERS];
    fp_t h2o_number_density[NUM_LAYERS];
    fp_t partial_pressure[NUM_LAYERS];
};


/*Initializes the atmospheric column test input data.*/
static struct Atmosphere initialize_atmosphere()
{
    struct Atmosphere atmosphere = {
        .level_pressure = {
            5.8692e-5, 0.002911421, 0.168689857, 0.486686405, 0.806365652,
            0.9994064
        },
        .level_xh2o = {
            2.04233e-6, 4.4842055e-06, 1.9884803e-05, 0.00075805407, 0.004059154,
            0.00771232
        },
        .pressure = {
            9.8692e-5, 0.00572415, 0.331655564, 0.641717246, 0.971014064
        },
        .temperature = {
            230.92, 236.24, 241.93, 275.15, 288.99
        },
        .xh2o = {
            4.072945e-06, 4.895466e-06, 3.487414e-05, 0.001481234, 0.006637074
        }
    };
    calc_number_densities(NUM_LAYERS, atmosphere.level_pressure, atmosphere.number_density);
    calc_partial_pressures_and_number_densities(NUM_LAYERS, atmosphere.level_pressure,
                                                atmosphere.level_xh2o, atmosphere.number_density,
                                                atmosphere.partial_pressure, atmosphere.h2o_number_density);
    return atmosphere;
}


/*Initializes a spectral grid object.*/
void initialize_spectral_bins(SpectralBins_t * bins)
{
    create_spectral_bins(bins, NUM_LAYERS, 475., NUM_GRID_POINTS, 1.,
                         2., HOST_ONLY);
}


/*Helper structure that mimics a set of molecular transitions.*/
struct Transitions
{
    fp_t center[NUM_TRANSITIONS];
    fp_t delta[NUM_TRANSITIONS];
    fp_t energy[NUM_TRANSITIONS];
    fp_t gamma_foreign[NUM_TRANSITIONS];
    fp_t gamma_self[NUM_TRANSITIONS];
    int isotopologue[NUM_TRANSITIONS];
    fp_t mass;
    int molecule[NUM_TRANSITIONS];
    fp_t n[NUM_TRANSITIONS];
    fp_t strength[NUM_TRANSITIONS];
};


/*Initializes  the molecular transitions test data.*/
static struct Transitions initialize_transitions()
{
    struct Transitions transitions = {
        .center = {
            500.018237, 500.03511, 500.436538, 500.621743, 500.659494
        },
        .delta = {
            0.00081, 0.0134, -0.01107, 0.01122, -0.0155
        },
        .energy = {
            4195.8181, 2248.0645, 4052.837, 1774.7511, 4006.0743
        },
        .gamma_foreign = {
            0.083, 0.0474, 0.0858, 0.0526, 0.0397
        },
        .gamma_self = {
            0.371, 0.268, 0.352, 0.284, 0.208
        },
        .isotopologue = {
            1, 1, 1, 1, 1
        },
        .mass = 18.010565/6.023e23,
        .molecule = {
            1, 1, 1, 1, 1
        },
        .n = {
            0.69, 0.49, 0.72, 0.53, 0.45
        },
        .strength = {
            1.356e-30, 2.363e-24, 1.309e-30, 3.416e-23, 2.394e-30
        }
    };
    return transitions;
}


int setup_optical_depth(fp_t * alpha, fp_t * gamma, fp_t * strength, fp_t * vnn)
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    rc_check(calc_line_centers(NUM_TRANSITIONS, NUM_LAYERS, transitions.center,
                               transitions.delta, atmosphere.pressure, vnn));
    rc_check(calc_doppler_hw(NUM_TRANSITIONS, NUM_LAYERS, transitions.mass,
                             vnn, atmosphere.temperature, alpha));
    rc_check(calc_lorentz_hw(NUM_TRANSITIONS, NUM_LAYERS, transitions.n,
                             transitions.gamma_foreign, transitions.gamma_self,
                             atmosphere.temperature, atmosphere.pressure,
                             atmosphere.partial_pressure, gamma));
    fp_t q[NUM_LAYERS*NUM_ISOTOPOLOGUES];
    rc_check(calc_partition_functions(NUM_LAYERS, transitions.molecule[0],
                                      NUM_ISOTOPOLOGUES, atmosphere.temperature, q));
    rc_check(calc_line_strengths(NUM_TRANSITIONS, NUM_LAYERS, NUM_ISOTOPOLOGUES,
                                 transitions.isotopologue, transitions.strength,
                                 transitions.center, transitions.energy,
                                 atmosphere.temperature, q, strength));
    return GRTCODE_SUCCESS;
}


/*Adjusts the line center at the input pressure.*/
int test_calc_line_centers()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(calc_line_centers(NUM_TRANSITIONS, NUM_LAYERS, transitions.center,
                               transitions.delta, atmosphere.pressure, vnn));
    fp_t vnn_ref[NUM_LAYERS*NUM_TRANSITIONS] = {
        500.01823708, 500.03511132, 500.43653691, 500.62174411, 500.65949247,
        500.01824164, 500.03518670, 500.43647463, 500.62180722, 500.65940528,
        500.01850564, 500.03955418, 500.43286657, 500.62546418, 500.65435334,
        500.01875679, 500.04370901, 500.42943419, 500.62894307, 500.64954738,
        500.01902352, 500.04812159, 500.42578887, 500.63263778, 500.64444328
    };
    rc_check(check_array(vnn, vnn_ref, NUM_LAYERS*NUM_TRANSITIONS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculates the partition function.*/
int test_calc_partition_functions()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    fp_t q[NUM_LAYERS*NUM_ISOTOPOLOGUES];
    rc_check(calc_partition_functions(NUM_LAYERS, transitions.molecule[0],
                                      NUM_ISOTOPOLOGUES, atmosphere.temperature, q));
    fp_t q_ref[NUM_LAYERS*NUM_ISOTOPOLOGUES] = {
        0.00829304, 0.00822385, 0.00137605, 0.00167721, 0.00165653, 0.00027761,
        0.00141451, 0.00139350, 0.00023390, 0.00801678, 0.00794978, 0.00133023,
        0.00162114, 0.00160113, 0.00026832, 0.00136700, 0.00134669, 0.00022604,
        0.00773788, 0.00767377, 0.00128397, 0.00156452, 0.00154519, 0.00025894,
        0.00131902, 0.00129942, 0.00021811, 0.00638825, 0.00633500, 0.00105998,
        0.00129043, 0.00127448, 0.00021353, 0.00108669, 0.00107048, 0.00017968,
        0.00593679, 0.00588716, 0.00098513, 0.00119881, 0.00118398, 0.00019835,
        0.00100896, 0.00099389, 0.00016683
    };
    rc_check(check_array(q, q_ref, NUM_LAYERS*NUM_ISOTOPOLOGUES, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculates the line strengths.*/
int test_calc_line_strengths()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    fp_t q[NUM_LAYERS*NUM_ISOTOPOLOGUES];
    rc_check(calc_partition_functions(NUM_LAYERS, transitions.molecule[0],
                                      NUM_ISOTOPOLOGUES, atmosphere.temperature, q));
    fp_t strength[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(calc_line_strengths(NUM_TRANSITIONS, NUM_LAYERS, NUM_ISOTOPOLOGUES,
                                 transitions.isotopologue, transitions.strength,
                                 transitions.center, transitions.energy,
                                 atmosphere.temperature, q, strength));
    fp_t strength_ref[NUM_LAYERS*NUM_TRANSITIONS] = {
        4.76159074e-44, 1.54670896e-32, 1.12041472e-43, 4.26855439e-30,
        2.74237586e-43, 8.26505087e-44, 2.04274624e-32, 1.90617316e-43,
        5.27531963e-30, 4.63513576e-43, 1.44965951e-43, 2.71062650e-32,
        3.27560268e-43, 6.54130472e-30, 7.91195631e-43, 2.37780419e-42,
        1.09803346e-31, 4.84878070e-42, 1.88643258e-29, 1.13253653e-41,
        6.25277511e-42, 1.77288209e-31, 1.23022854e-41, 2.70544514e-29,
        2.84003241e-41
    };
    rc_check(check_array(strength, strength_ref, NUM_LAYERS*NUM_TRANSITIONS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculates the lorentzian halfwidth at half-maximum.*/
int test_calc_lorentz_hw()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    fp_t gamma[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(calc_lorentz_hw(NUM_TRANSITIONS, NUM_LAYERS, transitions.n,
                             transitions.gamma_foreign, transitions.gamma_self,
                             atmosphere.temperature, atmosphere.pressure,
                             atmosphere.partial_pressure, gamma));
    fp_t gamma_ref[NUM_LAYERS*NUM_TRANSITIONS] = {
        9.72401695e-06, 5.28455223e-06, 1.01270053e-05, 5.92274072e-06,
        4.38225103e-06, 5.55517338e-04, 3.03335436e-04, 5.78109256e-04,
        3.39643874e-04, 2.51755074e-04, 3.16867913e-02, 1.73894432e-02,
        3.29492023e-02, 1.94512810e-02, 1.44448846e-02, 5.65140448e-02,
        3.19019347e-02, 5.84940512e-02, 3.54823934e-02, 2.66136084e-02,
        8.35112769e-02, 4.77690511e-02, 8.62183407e-02, 5.29874600e-02,
        3.98813419e-02
    };
    rc_check(check_array(gamma, gamma_ref, NUM_LAYERS*NUM_TRANSITIONS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculates the doppler halfwidth at half-maximum.*/
int test_calc_doppler_hw()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    struct Transitions transitions = initialize_transitions();
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(calc_line_centers(NUM_TRANSITIONS, NUM_LAYERS, transitions.center,
                               transitions.delta, atmosphere.pressure, vnn));
    fp_t alpha[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(calc_doppler_hw(NUM_TRANSITIONS, NUM_LAYERS, transitions.mass,
                             vnn, atmosphere.temperature, alpha));
    fp_t alpha_ref[NUM_LAYERS*NUM_TRANSITIONS] = {
        0.00064122, 0.00064125, 0.00064176, 0.00064200, 0.00064205, 0.00064857,
        0.00064859, 0.00064911, 0.00064935, 0.00064940, 0.00065633, 0.00065636,
        0.00065688, 0.00065713, 0.00065717, 0.00069994, 0.00069998, 0.00070052,
        0.00070080, 0.00070083, 0.00071733, 0.00071737, 0.00071792, 0.00071821,
        0.00071823
    };
    rc_check(check_array(alpha, alpha_ref, NUM_LAYERS*NUM_TRANSITIONS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


int test_sort_lines()
{
    struct Transitions transitions = initialize_transitions();
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t snn[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t gamma[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t alpha[NUM_LAYERS*NUM_TRANSITIONS];
    int i;
    for (i=0; i<NUM_LAYERS; ++i)
    {
        int j;
        for (j=0; j<NUM_TRANSITIONS; ++j)
        {
            vnn[i*NUM_TRANSITIONS + j] = transitions.center[j];
            snn[i*NUM_TRANSITIONS + j] = transitions.strength[j];
            gamma[i*NUM_TRANSITIONS + j] = transitions.gamma_self[j];
            alpha[i*NUM_TRANSITIONS + j] = transitions.gamma_foreign[j];
        }
        if (i < NUM_TRANSITIONS)
        {
            int swapping_index = i*NUM_TRANSITIONS + NUM_TRANSITIONS - 1 - i;
            fp_t x = vnn[i*NUM_TRANSITIONS];
            vnn[i*NUM_TRANSITIONS] = vnn[swapping_index];
            vnn[swapping_index] = x;
            x = snn[i*NUM_TRANSITIONS];
            snn[i*NUM_TRANSITIONS] = snn[swapping_index];
            snn[swapping_index] = x;
            x = gamma[i*NUM_TRANSITIONS];
            gamma[i*NUM_TRANSITIONS] = gamma[swapping_index];
            gamma[swapping_index] = x;
            x = alpha[i*NUM_TRANSITIONS];
            alpha[i*NUM_TRANSITIONS] = alpha[swapping_index];
            alpha[swapping_index] = x;
        }
    }
    rc_check(sort_lines(NUM_TRANSITIONS, NUM_LAYERS, vnn, snn, gamma, alpha));
    for (i=0; i<NUM_LAYERS; ++i)
    {
        int j;
        for (j=0; j<NUM_TRANSITIONS; ++j)
        {
            if (vnn[i*NUM_TRANSITIONS + j] != transitions.center[j] ||
                snn[i*NUM_TRANSITIONS + j] != transitions.strength[j] ||
                gamma[i*NUM_TRANSITIONS + j] != transitions.gamma_self[j] ||
                alpha[i*NUM_TRANSITIONS + j] != transitions.gamma_foreign[j])
            {
                return GRTCODE_VALUE_ERR;
            }
        }
    }
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_bin_sweep()
{
    return GRTCODE_SUCCESS;
    struct Transitions transitions = initialize_transitions();
    fp_t alpha[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t gamma[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t strength[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(setup_optical_depth(alpha, gamma, strength, vnn));
    SpectralBins_t bins;
    initialize_spectral_bins(&bins);
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_optical_depth_bin_sweep(NUM_TRANSITIONS, NUM_LAYERS, vnn, strength, gamma,
                                          alpha, transitions.n, bins, tau));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS] = {
    };
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_line_sweep()
{
    return GRTCODE_SUCCESS;
    struct Transitions transitions = initialize_transitions();
    fp_t alpha[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t gamma[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t strength[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(setup_optical_depth(alpha, gamma, strength, vnn));
    SpectralBins_t bins;
    initialize_spectral_bins(&bins);
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_optical_depth_line_sweep(NUM_TRANSITIONS, NUM_LAYERS, vnn, strength, gamma,
                                           alpha, transitions.n, bins, tau));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS] = {
    };
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int test_calc_optical_depth_line_sample()
{
    return GRTCODE_SUCCESS;
    struct Transitions transitions = initialize_transitions();
    fp_t alpha[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t gamma[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t strength[NUM_LAYERS*NUM_TRANSITIONS];
    fp_t vnn[NUM_LAYERS*NUM_TRANSITIONS];
    rc_check(setup_optical_depth(alpha, gamma, strength, vnn));
    SpectralBins_t bins;
    initialize_spectral_bins(&bins);
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_optical_depth_line_sample(NUM_TRANSITIONS, NUM_LAYERS, vnn, strength, gamma,
                                            alpha, transitions.n, bins, tau, NULL, NULL));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS] = {
    };
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int test_calc_water_vapor_ctm_optical_depth()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    fp_t T0F[NUM_GRID_POINTS] = {
        3.004e-03, 3.020e-03, 3.030e-03, 3.111e-03, 3.130e-03,
        3.161e-03, 3.180e-03, 3.237e-03, 2.463e-03, 2.473e-03,
        2.484e-03, 2.515e-03, 2.477e-03, 2.490e-03, 2.502e-03,
        2.518e-03, 2.530e-03, 2.711e-03, 2.754e-03, 2.744e-03,
        2.898e-03, 2.910e-03, 2.928e-03, 2.437e-03, 2.446e-03,
        2.459e-03, 2.697e-03, 2.697e-03, 2.712e-03, 2.728e-03,
        2.745e-03, 2.762e-03, 2.769e-03, 2.789e-03, 2.661e-03,
        2.681e-03, 2.694e-03, 2.706e-03, 2.721e-03, 2.736e-03,
        2.749e-03, 2.766e-03, 2.771e-03, 2.658e-03, 2.675e-03,
        2.686e-03, 2.737e-03, 2.805e-03, 2.856e-03, 2.873e-03
    };
    fp_t CF[NUM_GRID_POINTS] = {
        2.32e-23, 2.29e-23, 2.27e-23, 2.22e-23, 2.19e-23,
        2.16e-23, 2.13e-23, 2.09e-23, 2.25e-23, 2.22e-23,
        2.19e-23, 2.16e-23, 2.14e-23, 2.10e-23, 2.07e-23,
        2.04e-23, 2.02e-23, 1.96e-23, 1.93e-23, 1.91e-23,
        1.86e-23, 1.83e-23, 1.81e-23, 1.90e-23, 1.87e-23,
        1.84e-23, 1.78e-23, 1.76e-23, 1.73e-23, 1.71e-23,
        1.68e-23, 1.66e-23, 1.64e-23, 1.61e-23, 1.63e-23,
        1.61e-23, 1.59e-23, 1.56e-23, 1.54e-23, 1.52e-23,
        1.50e-23, 1.48e-23, 1.48e-23, 1.46e-23, 1.44e-23,
        1.42e-23, 1.40e-23, 1.37e-23, 1.35e-23, 1.34e-23
    };
    fp_t T0S[NUM_GRID_POINTS] = {
        1.459e-02, 1.462e-02, 1.465e-02, 1.468e-02, 1.472e-02,
        1.475e-02, 1.478e-02, 1.482e-02, 1.483e-02, 1.486e-02,
        1.490e-02, 1.493e-02, 1.497e-02, 1.500e-02, 1.504e-02,
        1.508e-02, 1.511e-02, 1.515e-02, 1.518e-02, 1.521e-02,
        1.525e-02, 1.529e-02, 1.532e-02, 1.534e-02, 1.538e-02,
        1.542e-02, 1.545e-02, 1.548e-02, 1.551e-02, 1.555e-02,
        1.558e-02, 1.561e-02, 1.565e-02, 1.568e-02, 1.571e-02,
        1.575e-02, 1.578e-02, 1.580e-02, 1.583e-02, 1.585e-02,
        1.588e-02, 1.591e-02, 1.594e-02, 1.596e-02, 1.599e-02,
        1.602e-02, 1.605e-02, 1.608e-02, 1.611e-02, 1.614e-02
    };
    fp_t CS[NUM_GRID_POINTS] = {
        4.96e-21, 4.91e-21, 4.86e-21, 4.82e-21, 4.77e-21,
        4.72e-21, 4.68e-21, 4.64e-21, 4.60e-21, 4.56e-21,
        4.52e-21, 4.48e-21, 4.43e-21, 4.39e-21, 4.35e-21,
        4.30e-21, 4.27e-21, 4.22e-21, 4.18e-21, 4.15e-21,
        4.11e-21, 4.07e-21, 4.03e-21, 3.99e-21, 3.95e-21,
        3.91e-21, 3.88e-21, 3.84e-21, 3.81e-21, 3.77e-21,
        3.74e-21, 3.70e-21, 3.67e-21, 3.63e-21, 3.60e-21,
        3.57e-21, 3.53e-21, 3.50e-21, 3.47e-21, 3.44e-21,
        3.41e-21, 3.38e-21, 3.35e-21, 3.32e-21, 3.28e-21,
        3.25e-21, 3.22e-21, 3.19e-21, 3.16e-21, 3.13e-21
    };
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    memset(tau, 0, sizeof(*tau)*NUM_LAYERS*NUM_GRID_POINTS);
    rc_check(calc_water_vapor_ctm_optical_depth(NUM_GRID_POINTS, NUM_LAYERS,
                                                tau, CS, atmosphere.temperature,
                                                atmosphere.partial_pressure,
                                                atmosphere.h2o_number_density, T0S, CF,
                                                atmosphere.pressure, T0F));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS] = {
        7.31326475e-10, 7.22688704e-10, 7.16840542e-10, 7.04941728e-10,
        6.96354323e-10, 6.88256871e-10, 6.79642615e-10, 6.69513978e-10,
        6.84825627e-10, 6.76234144e-10, 6.67683990e-10, 6.59947042e-10,
        6.52272875e-10, 6.40798387e-10, 6.32250840e-10, 6.23810165e-10,
        6.18240699e-10, 6.07129177e-10, 5.99590455e-10, 5.93081123e-10,
        5.83517617e-10, 5.74681174e-10, 5.69095478e-10, 5.78220509e-10,
        5.69541123e-10, 5.60992976e-10, 5.51345736e-10, 5.45193037e-10,
        5.36572736e-10, 5.30960105e-10, 5.22373023e-10, 5.16755047e-10,
        5.10848059e-10, 5.02270209e-10, 5.04146202e-10, 4.98681791e-10,
        4.92935855e-10, 4.84162389e-10, 4.78493074e-10, 4.72803268e-10,
        4.67051430e-10, 4.61405925e-10, 4.61459917e-10, 4.52060448e-10,
        4.46393881e-10, 4.40586342e-10, 4.35870073e-10, 4.28545327e-10,
        4.23741823e-10, 4.21052925e-10, 9.44781330e-06, 9.33755183e-06,
        9.26132011e-06, 9.10991504e-06, 9.00049569e-06, 8.89614914e-06,
        8.78729154e-06, 8.65885920e-06, 8.87274882e-06, 8.76412594e-06,
        8.65627742e-06, 8.55714356e-06, 8.46070860e-06, 8.31693171e-06,
        8.20903157e-06, 8.10103080e-06, 8.03018696e-06, 7.88346804e-06,
        7.78639574e-06, 7.70504075e-06, 7.57968828e-06, 7.46842183e-06,
        7.39600846e-06, 7.51991373e-06, 7.41035002e-06, 7.30222785e-06,
        7.17394368e-06, 7.09514314e-06, 6.98716886e-06, 6.91463494e-06,
        6.80700111e-06, 6.73408722e-06, 6.65953262e-06, 6.55064301e-06,
        6.57438826e-06, 6.50471577e-06, 6.43029665e-06, 6.32022319e-06,
        6.24795486e-06, 6.17518019e-06, 6.10196943e-06, 6.02992923e-06,
        6.02730679e-06, 5.91319522e-06, 5.83954725e-06, 5.76565423e-06,
        5.70389632e-06, 5.61009037e-06, 5.54732632e-06, 5.51108588e-06,
        3.46246040e-02, 3.42252490e-02, 3.39427275e-02, 3.33977856e-02,
        3.30020750e-02, 3.26206042e-02, 3.22307013e-02, 3.17695498e-02,
        3.25996597e-02, 3.22103974e-02, 3.18246903e-02, 3.14646149e-02,
        3.11202006e-02, 3.06103126e-02, 3.02241401e-02, 2.98319805e-02,
        2.95763144e-02, 2.90289647e-02, 2.86747116e-02, 2.83864283e-02,
        2.79220484e-02, 2.75249241e-02, 2.72583738e-02, 2.77270827e-02,
        2.73350357e-02, 2.69473628e-02, 2.64665737e-02, 2.61800270e-02,
        2.57972316e-02, 2.55309599e-02, 2.51491050e-02, 2.48804768e-02,
        2.46135135e-02, 2.42218552e-02, 2.43046346e-02, 2.40525672e-02,
        2.37789905e-02, 2.33882571e-02, 2.31270183e-02, 2.28631931e-02,
        2.25989572e-02, 2.23382846e-02, 2.23155476e-02, 2.19240100e-02,
        2.16523933e-02, 2.13858190e-02, 2.11566438e-02, 2.08166210e-02,
        2.05842143e-02, 2.04455085e-02, 4.78057745e-01, 4.72651417e-01,
        4.68419003e-01, 4.61339899e-01, 4.55970834e-01, 4.50622282e-01,
        4.45600181e-01, 4.39596230e-01, 4.53143960e-01, 4.48121961e-01,
        4.43139654e-01, 4.38223912e-01, 4.33770200e-01, 4.27552229e-01,
        4.22560643e-01, 4.17184308e-01, 4.13749303e-01, 4.05572118e-01,
        4.00665117e-01, 3.97104283e-01, 3.90326717e-01, 3.85280653e-01,
        3.81434054e-01, 3.88441262e-01, 3.83409764e-01, 3.78392738e-01,
        3.71149135e-01, 3.67221207e-01, 3.62556322e-01, 3.58724533e-01,
        3.54060513e-01, 3.50193465e-01, 3.46718077e-01, 3.41649435e-01,
        3.42412766e-01, 3.38988743e-01, 3.35096110e-01, 3.30368696e-01,
        3.26884704e-01, 3.23368532e-01, 3.19871307e-01, 3.16386202e-01,
        3.15266819e-01, 3.11271596e-01, 3.07376532e-01, 3.03865159e-01,
        3.00491918e-01, 2.95963965e-01, 2.92575795e-01, 2.90276198e-01,
        1.30972953e+00, 1.29511346e+00, 1.28293466e+00, 1.26511623e+00,
        1.25054406e+00, 1.23596137e+00, 1.22281775e+00, 1.20738721e+00,
        1.23797384e+00, 1.22482753e+00, 1.21172873e+00, 1.19865305e+00,
        1.18632941e+00, 1.17072420e+00, 1.15761303e+00, 1.14301866e+00,
        1.13380307e+00, 1.11238501e+00, 1.09930579e+00, 1.09000056e+00,
        1.07238820e+00, 1.05922919e+00, 1.04850574e+00, 1.06317333e+00,
        1.05002051e+00, 1.03687613e+00, 1.01850561e+00, 1.00772069e+00,
        9.96015547e-01, 9.85312206e-01, 9.73605618e-01, 9.62857969e-01,
        9.53615505e-01, 9.40407863e-01, 9.40587113e-01, 9.31376591e-01,
        9.20605742e-01, 9.08832491e-01, 8.99561287e-01, 8.90250990e-01,
        8.80969117e-01, 8.71694691e-01, 8.67308969e-01, 8.57667179e-01,
        8.46887115e-01, 8.37591271e-01, 8.28389368e-01, 8.16761024e-01,
        8.07546968e-01, 8.00712224e-01
    };
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


int test_calc_ozone_ctm_optical_depth()
{
    return GRTCODE_SUCCESS;
    struct Atmosphere atmosphere = initialize_atmosphere();
    fp_t cross_section[NUM_GRID_POINTS] = {
        9.666e-25, 1.018e-24, 1.069e-24, 1.120e-24, 1.171e-24,
        1.222e-24, 1.273e-24, 1.248e-24, 1.223e-24, 1.198e-24,
        1.173e-24, 1.148e-24, 1.123e-24, 1.098e-24, 1.173e-24,
        1.247e-24, 1.322e-24, 1.396e-24, 1.471e-24, 1.545e-24,
        1.620e-24, 1.588e-24, 1.556e-24, 1.524e-24, 1.492e-24,
        1.460e-24, 1.428e-24, 1.396e-24, 1.363e-24, 1.329e-24,
        1.296e-24, 1.263e-24, 1.230e-24, 1.196e-24, 1.163e-24,
        1.183e-24, 1.203e-24, 1.223e-24, 1.242e-24, 1.262e-24,
        1.282e-24, 1.302e-24, 1.305e-24, 1.308e-24, 1.311e-24,
        1.314e-24, 1.317e-24, 1.320e-24, 1.323e-24, 1.326e-24
    };
    fp_t n[NUM_LAYERS];
    rc_check(calc_number_densities(NUM_LAYERS, atmosphere.level_pressure, n));
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_ozone_ctm_optical_depth(NUM_GRID_POINTS, NUM_LAYERS, cross_section,
                                          n, tau));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


int test_interpolate()
{
    return GRTCODE_SUCCESS;
    SpectralBins_t bins;
    initialize_spectral_bins(&bins);
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(interpolate(bins, tau));
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int test_interpolate_last_bin()
{
    return GRTCODE_SUCCESS;
    SpectralBins_t bins;
    initialize_spectral_bins(&bins);
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(interpolate_last_bin(bins, tau));
    rc_check(destroy_spectral_bins(&bins));
    return GRTCODE_SUCCESS;
}


int test_calc_cfc_optical_depth()
{
    return GRTCODE_SUCCESS;
    struct Atmosphere atmosphere = initialize_atmosphere();
    fp_t cross_section[NUM_GRID_POINTS] = {
        3.526e-20, 4.106e-20, 3.763e-20, 3.038e-20, 2.733e-20,
        2.977e-20, 3.411e-20, 3.590e-20, 3.492e-20, 3.629e-20,
        3.916e-20, 3.599e-20, 3.319e-20, 4.220e-20, 5.096e-20,
        4.686e-20, 4.331e-20, 4.773e-20, 4.750e-20, 4.254e-20,
        4.366e-20, 4.605e-20, 4.383e-20, 4.443e-20, 4.539e-20,
        3.914e-20, 3.806e-20, 4.557e-20, 4.192e-20, 3.129e-20,
        3.384e-20, 3.873e-20, 3.127e-20, 2.546e-20, 3.296e-20,
        4.108e-20, 3.645e-20, 3.066e-20, 3.830e-20, 4.621e-20,
        4.389e-20, 4.293e-20, 4.337e-20, 3.951e-20, 3.949e-20,
        4.354e-20, 4.496e-20, 4.414e-20, 4.310e-20, 4.438e-20
    };
    fp_t n[NUM_LAYERS];
    rc_check(calc_number_densities(NUM_LAYERS, atmosphere.level_pressure, n));
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_cfc_optical_depth(NUM_GRID_POINTS, NUM_LAYERS, n, atmosphere.level_xh2o,
                                    cross_section, tau));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


int test_calc_cia_optical_depth()
{
    struct Atmosphere atmosphere = initialize_atmosphere();
    fp_t cross_section[NUM_GRID_POINTS] = {
        2.163e-47, 2.118e-47, 2.072e-47, 2.027e-47, 1.979e-47,
        1.934e-47, 1.892e-47, 1.853e-47, 1.814e-47, 1.775e-47,
        1.733e-47, 1.691e-47, 1.650e-47, 1.611e-47, 1.572e-47,
        1.535e-47, 1.500e-47, 1.466e-47, 1.431e-47, 1.397e-47,
        1.364e-47, 1.332e-47, 1.300e-47, 1.269e-47, 1.239e-47,
        1.209e-47, 1.180e-47, 1.151e-47, 1.121e-47, 1.094e-47,
        1.067e-47, 1.040e-47, 1.016e-47, 9.906e-48, 9.662e-48,
        9.423e-48, 9.190e-48, 8.964e-48, 8.742e-48, 8.518e-48,
        8.308e-48, 8.102e-48, 7.894e-48, 7.700e-48, 7.509e-48,
        7.324e-48, 7.136e-48, 6.960e-48, 6.781e-48, 6.614e-48
    };
    fp_t tau[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(calc_cia_optical_depth(NUM_GRID_POINTS, NUM_LAYERS, atmosphere.level_pressure,
                                    atmosphere.temperature, atmosphere.level_xh2o,
                                    atmosphere.level_xh2o, cross_section, tau));
    fp_t tau_ref[NUM_LAYERS*NUM_GRID_POINTS];
    rc_check(check_array(tau, tau_ref, NUM_LAYERS*NUM_GRID_POINTS, 1.e-8, 0.));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_calc_line_centers,
            "test_calc_line_centers",
            GRTCODE_SUCCESS
        },
        {
            test_calc_partition_functions,
            "test_calc_partition_functions",
            GRTCODE_SUCCESS
        },
        {
            test_calc_line_strengths,
            "test_calc_line_strengths",
            GRTCODE_SUCCESS
        },
        {
            test_calc_lorentz_hw,
            "test_calc_lorentz_hw",
            GRTCODE_SUCCESS
        },
        {
            test_calc_doppler_hw,
            "test_calc_doppler_hw",
            GRTCODE_SUCCESS
        },
        {
            test_sort_lines,
            "test_sort_lines",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_bin_sweep,
            "test_calc_optical_depth_bin_sweep",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_line_sweep,
            "test_calc_optical_depth_line_sweep",
            GRTCODE_SUCCESS
        },
        {
            test_calc_optical_depth_line_sample,
            "test_calc_optical_depth_line_sample",
            GRTCODE_SUCCESS
        },
        {
            test_calc_water_vapor_ctm_optical_depth,
            "test_calc_water_vapor_ctm_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_calc_ozone_ctm_optical_depth,
            "test_calc_ozone_ctm_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate,
            "test_interpolate",
            GRTCODE_SUCCESS
        },
        {
            test_interpolate_last_bin,
            "test_interpolate_last_bin",
            GRTCODE_SUCCESS
        },
        {
            test_calc_cfc_optical_depth,
            "test_calc_cfc_optical_depth",
            GRTCODE_SUCCESS
        },
        {
            test_calc_cia_optical_depth,
            "test_calc_cia_optical_depth",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_kernels", NUM_TESTS, tests);
}
