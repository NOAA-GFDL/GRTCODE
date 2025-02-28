#include "grtcode_utilities.h"
#include "molecules.h"
#include "test_harness.h"


#define NUM_TESTS 3


int test_create_molecule()
{
    Molecule_t molecule;
    int id = H2O;
    char * hitran_path = "";
    double w0 = 1.;
    double wn = 1000.;
    int num_layers = 10;
    Device_t device;
    rc_check(create_molecule(&molecule, id, hitran_path, w0, wn, num_layers, device));
    return GRTCODE_SUCCESS;
}


int test_free_molecule()
{
    Molecule_t molecule;
    rc_check(free_molecule(&molecule));
    return GRTCODE_SUCCESS;
}


int test_molecule_hash()
{
    int id = H2O;
    int hash;
    rc_check(molecule_hash(id, &hash));
    return GRTCODE_SUCCESS;
}


int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_create_molecule,
            "test_create_molecule",
            GRTCODE_SUCCESS
        },
        {
            test_free_molecule,
            "test_free_molecule",
            GRTCODE_SUCCESS
        },
        {
            test_molecule_hash,
            "test_molecule_hash",
            GRTCODE_SUCCESS
        }
    };
    return test_harness("test_molecule", NUM_TESTS, tests);
}
