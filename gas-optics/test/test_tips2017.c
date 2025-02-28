#include "grtcode_utilities.h"
#include "molecules.h"
#include "test_harness.h"
#include "tips2017.h"


#define NUM_TESTS 6
#define TEMPERATURE 275.234324
#define ISOTOPOLOGUE_VALUE 1


/*Helper function that runs the test for an input molecule.*/
int helper(int molecule, fp_t result)
{
    rc_check(inittips_d());
    fp_t qt = Q(molecule, TEMPERATURE, ISOTOPOLOGUE_VALUE);
    if (check_floating_point(qt, result, 1.e-9, 0.))
    {
        return GRTCODE_SUCCESS;
    }
    return GRTCODE_VALUE_ERR;
}


/*Run the initialization function.*/
int test_inittips_d()
{
    rc_check(inittips_d());
    return GRTCODE_SUCCESS;
}


/*Calculate the water total internal partition function value.*/
int test_Q_H2O()
{
    return helper(H2O, 156.6091754);
}


/*Calculate the carbon dioxide total internal partition function value.*/
int test_Q_CO2()
{
    return helper(CO2, 261.25798746);
}


/*Calculate the methane total internal partition function value.*/
int test_Q_CH4()
{
    return helper(CH4, 528.2642260800001);
}


/*Calculate the nitrous oxide total internal partition function value.*/
int test_Q_N2O()
{
    return helper(N2O, 4524.7762498);
}


/*Calculate the ozone total internal partition function value.*/
int test_Q_O3()
{
    return helper(O3, 3087.3115616000005);
}


/*Runs the tests.*/
int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_inittips_d,
            "test_inittips_d",
            GRTCODE_SUCCESS
        },
        {
            test_Q_H2O,
            "test_Q_H2O",
            GRTCODE_SUCCESS
        },
        {
            test_Q_CO2,
            "test_Q_CO2",
            GRTCODE_SUCCESS
        },
        {
            test_Q_CH4,
            "test_Q_CH4",
            GRTCODE_SUCCESS
        },
        {
            test_Q_N2O,
            "test_Q_N2O",
            GRTCODE_SUCCESS
        },
        {
            test_Q_O3,
            "test_Q_O3",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_tips2017", NUM_TESTS, tests);
}
