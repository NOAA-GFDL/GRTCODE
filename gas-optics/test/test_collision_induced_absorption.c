#include "collision_induced_absorption.h"
#include "collision_induced_absorption-internal.h"
#include "grtcode_utilities.h"
#include "test_harness.h"


#define NUM_TESTS 2


int test_get_collision_induced_cross_sections()
{
    Device_t device;
    SpectralGrid_t grid;
    int id[2] = {};
    char * filepath = "";
    CollisionInducedAbsorption_t cia;
    rc_check(get_collision_induced_cross_sections(&cia, id, filepath, grid, device));
    return GRTCODE_SUCCESS;
}


int test_free_collision_induced_cross_sections()
{
    CollisionInducedAbsorption_t cia;
    rc_check(free_collision_induced_cross_sections(&cia));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_get_collision_induced_cross_sections,
            "test_get_collision_induced_cross_sections",
            GRTCODE_SUCCESS
        },
        {
            test_free_collision_induced_cross_sections,
            "test_free_collision_induced_cross_sections",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_collision_induced_absorption", NUM_TESTS, tests);
}
