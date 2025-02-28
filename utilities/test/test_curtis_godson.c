#include "circ1.h"
#include "curtis_godson.h"
#include "floating_point_type.h"
#include "return_codes.h"
#include "test_harness.h"


#define NUM_TESTS 3


/*Calculate integrated number densities.*/
int test_calc_number_densities()
{
    fp_t p[NUM_LEVELS];
    int i;
    for (i=0; i<num_levels; ++i)
    {
        fp_t const mbtoatm = 0.000986923;
        p[i] = level_pressure[i]*mbtoatm;
    }
    fp_t n[NUM_LAYERS];
    rc_check(calc_number_densities(num_layers, p, n));

    fp_t n_ref[NUM_LAYERS] = {
        1.27184116e+21, 1.90776174e+21, 1.48381468e+21, 1.69578821e+21,
        2.33170879e+21, 2.96762937e+21, 3.81552347e+21, 4.66341758e+21,
        6.14723226e+21, 8.05499399e+21, 1.08106498e+22, 1.42022263e+22,
        2.09853791e+22, 2.52248496e+22, 3.13720819e+22, 4.51503611e+22,
        7.01632372e+22, 6.46519255e+22, 8.01259929e+22, 1.01959266e+23,
        1.29515824e+23, 1.63855536e+23, 2.08581950e+23, 2.64966908e+23,
        3.75193141e+23, 4.74820699e+23, 6.06668232e+23, 4.87751084e+23,
        5.61941818e+23, 6.37192420e+23, 7.19650121e+23, 8.16522023e+23,
        4.48959928e+23, 4.75668593e+23, 5.03437125e+23, 5.29297895e+23,
        5.60881950e+23, 5.93525873e+23, 6.25745849e+23, 6.59873587e+23,
        6.93577378e+23, 7.31732612e+23, 7.71583635e+23, 8.13978340e+23,
        8.55313178e+23, 8.96012095e+23, 9.50277318e+23, 1.01111372e+24,
        1.08297275e+24, 1.13956968e+24, 5.86742720e+23, 6.04548497e+23,
        3.06725692e+23, 1.33543321e+23
    };
    rc_check(check_array(n, n_ref, num_layers, 1.e-9, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculate layer pressures and temperatures.*/
int test_calc_pressures_and_temperatures()
{
    fp_t p[NUM_LEVELS];
    int i;
    for (i=0; i<num_levels; ++i)
    {
        fp_t const mbtoatm = 0.000986923;
        p[i] = level_pressure[i]*mbtoatm;
    }
    fp_t pavg[NUM_LAYERS];
    fp_t tavg[NUM_LAYERS];
    rc_check(calc_pressures_and_temperatures(num_layers, p, level_temperature,
                                             pavg, tavg));

    fp_t pavg_ref[NUM_LAYERS] = {
        9.86923000e-05, 1.72711525e-04, 2.51665365e-04, 3.25684590e-04,
        4.19442275e-04, 5.42807650e-04, 7.00715330e-04, 8.98099930e-04,
        1.14976529e-03, 1.48038450e-03, 1.91956524e-03, 2.50184981e-03,
        3.32099590e-03, 4.39674197e-03, 5.71428417e-03, 7.49568019e-03,
        1.01801107e-02, 1.33185259e-02, 1.66888679e-02, 2.09277022e-02,
        2.63163018e-02, 3.31458090e-02, 4.18159275e-02, 5.28398574e-02,
        6.77423947e-02, 8.75302009e-02, 1.12706607e-01, 1.38184024e-01,
        1.62620237e-01, 1.90535354e-01, 2.22121825e-01, 2.57882980e-01,
        2.87342631e-01, 3.08867422e-01, 3.31660409e-01, 3.55701853e-01,
        3.81080578e-01, 4.07954491e-01, 4.36338397e-01, 4.66266837e-01,
        4.97774354e-01, 5.30954705e-01, 5.65950994e-01, 6.02861915e-01,
        6.41722008e-01, 6.82491797e-01, 7.25472293e-01, 7.71132286e-01,
        8.19881348e-01, 8.71620786e-01, 9.11808290e-01, 9.39540827e-01,
        9.60754737e-01, 9.71003932e-01
    };
    rc_check(check_array(pavg, pavg_ref, num_layers, 1.e-9, 0.));
    fp_t tavg_ref[NUM_LAYERS] = {
        230.480, 241.520, 249.865, 255.475, 260.840, 265.700, 269.320, 270.875, 269.095,
        264.305, 258.705, 253.170, 247.385, 241.685, 236.370, 231.430, 227.810, 225.720,
        224.230, 222.805, 221.345, 219.850, 218.350, 217.120, 210.935, 204.420, 204.110,
        206.635, 212.275, 219.350, 225.785, 230.955, 234.950, 238.320, 241.840, 245.770,
        249.355, 252.405, 255.600, 259.230, 262.780, 266.050, 268.900, 271.535, 275.035,
        278.440, 279.865, 278.975, 277.600, 279.725, 283.385, 285.885, 287.655, 288.625
    };
    rc_check(check_array(tavg, tavg_ref, num_layers, 1.e-9, 0.));
    return GRTCODE_SUCCESS;
}


/*Calculate partial pressures and number densities.*/
int test_calc_partial_pressures_and_number_densities()
{
    fp_t p[NUM_LEVELS];
    int i;
    for (i=0; i<num_levels; ++i)
    {
        fp_t const mbtoatm = 0.000986923;
        p[i] = level_pressure[i]*mbtoatm;
    }
    fp_t n[NUM_LAYERS];
    rc_check(calc_number_densities(num_layers, p, n));
    fp_t abundance[NUM_LEVELS];
    abundance[0] = H2O_abundance[0];
    abundance[num_layers] = H2O_abundance[num_layers - 1];
    for (i=1; i<num_layers; ++i)
    {
        abundance[i] = H2O_abundance[i - 1] + (H2O_abundance[i] - H2O_abundance[i - 1])*
                       (level_pressure[i] - layer_pressure[i - 1])/
                       (layer_pressure[i] - layer_pressure[i - 1]);
    }
    fp_t ps[NUM_LAYERS];
    fp_t ns[NUM_LAYERS];
    rc_check(calc_partial_pressures_and_number_densities(num_layers, p,
                                                         abundance, n,
                                                         ps, ns));

    fp_t ps_ref[NUM_LAYERS] = {
        4.13101186e-10, 7.81155351e-10, 1.21390593e-09, 1.61903859e-09,
        2.13002180e-09, 2.79873183e-09, 3.64671512e-09, 4.69836250e-09,
        6.02161073e-09, 7.71943518e-09, 9.90199574e-09, 1.27105251e-08,
        1.66144587e-08, 2.17351782e-08, 2.79739064e-08, 3.62756239e-08,
        4.84818845e-08, 6.23207525e-08, 7.67308646e-08, 9.42399415e-08,
        1.15182679e-07, 1.39887809e-07, 1.69947001e-07, 2.05324792e-07,
        2.58134561e-07, 3.27516798e-07, 4.03998659e-07, 5.17253656e-07,
        7.29885137e-07, 1.28857623e-06, 3.07516894e-06, 7.49570258e-06,
        1.12624880e-05, 1.14802593e-05, 1.43313841e-05, 2.18319289e-05,
        2.66083292e-05, 3.03553687e-05, 4.02441332e-05, 5.80282903e-05,
        9.47062574e-05, 1.56477200e-04, 2.96379097e-04, 5.79586095e-04,
        8.57072705e-04, 1.31046523e-03, 2.63019586e-03, 4.32640968e-03,
        5.28760622e-03, 5.75257097e-03, 6.07213400e-03, 6.24817060e-03,
        6.38677083e-03, 6.44850971e-03
    };
    rc_check(check_array(ps, ps_ref, num_layers, 1.e-9, 0.));
    fp_t ns_ref[NUM_LAYERS] = {
        5.31056512e+15, 8.59172063e+15, 7.15103954e+15, 8.42602533e+15,
        1.18360914e+16, 1.52973554e+16, 1.98540720e+16, 2.43946922e+16,
        3.21953490e+16, 4.20091229e+16, 5.57833273e+16, 7.21797768e+16,
        1.05023844e+17, 1.24725410e+17, 1.53612086e+17, 2.18573195e+17,
        3.34314246e+17, 3.02618443e+17, 3.68523235e+17, 4.59335492e+17,
        5.67235000e+17, 6.92052290e+17, 8.48307917e+17, 1.03067321e+18,
        1.42931213e+18, 1.77840803e+18, 2.17640371e+18, 1.82281653e+18,
        2.51436471e+18, 4.28072964e+18, 9.86791098e+18, 2.35357599e+19,
        1.75943031e+19, 1.76943392e+19, 2.17040864e+19, 3.24272717e+19,
        3.91700874e+19, 4.41255264e+19, 5.76293041e+19, 8.19779279e+19,
        1.31618533e+20, 2.15187942e+20, 4.02686073e+20, 7.80268524e+20,
        1.14138442e+21, 1.71613500e+21, 3.43339232e+21, 5.66495109e+21,
        6.98344328e+21, 7.52019020e+21, 3.90740597e+21, 4.02039167e+21,
        2.03901008e+21, 8.86871754e+20
    };
    rc_check(check_array(ns, ns_ref, num_layers, 1.e-9, 0.));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_calc_number_densities,
            "test_calc_number_densities",
            GRTCODE_SUCCESS
        },
        {
            test_calc_pressures_and_temperatures,
            "test_calc_pressures_and_temperatures",
            GRTCODE_SUCCESS
        },
        {
            test_calc_partial_pressures_and_number_densities,
            "test_calc_partial_pressures_and_number_densities",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_curtis_godson", NUM_TESTS, tests);
}
