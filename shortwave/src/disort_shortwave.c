#include <stdint.h>
#ifdef USE_DISORT
#include "cdisort.h"
#endif
#include "debug.h"
#include "disort_shortwave.h"
#include "grtcode_utilities.h"


/*Calculate upward and downward shortwave fluxes using DISORT.*/
EXTERN int disort_shortwave(Optics_t * const optics, fp_t const zen_dir,
                            fp_t * const surface_albedo, fp_t const total_solar_irradiance,
                            fp_t * const solar_flux, fp_t * const flux_up,
                            fp_t * const flux_down)
{
#ifdef USE_DISORT
    /*Configure DISORT.*/
    disort_state ds;
    ds.nlyr = optics->num_layers;
    ds.nstr = 16;
    ds.nmom = ds.nstr;
    ds.nphi = 0;
    ds.flag.planck = FALSE;
    ds.flag.onlyfl = TRUE;
    ds.accur = 0.;
    ds.flag.prnt[0] = FALSE;
    ds.flag.prnt[1] = FALSE;
    ds.flag.prnt[2] = FALSE;
    ds.flag.prnt[3] = FALSE;
    ds.flag.prnt[4] = FALSE;
    ds.flag.ibcnd = GENERAL_BC;
    ds.flag.usrang = FALSE;
    ds.flag.usrtau = FALSE;
    ds.flag.lamber = TRUE;

    /*Allocate memory.*/
    c_disort_state_alloc(&ds);
    disort_output out;
    c_disort_out_alloc(&ds, &out);
    int num_levels = optics->num_layers + 1;
    uint64_t n;
    for (n=0; n<(optics->grid.n); ++n)
    {
        int i;
        for (i=0; i<(optics->num_layers); ++i)
        {
            uint64_t o = i*(optics->grid.n) + n;
            ds.dtauc[i] = optics->tau[o];
            ds.ssalb[i] = optics->omega[o];
            c_getmom(HENYEY_GREENSTEIN, optics->g[o], ds.nmom,
                     &(ds.pmom[i*(ds.nmom_nstr+1)]));
        }
        ds.bc.umu0 = zen_dir;
        ds.bc.fbeam = solar_flux[n];
        ds.bc.phi0 = 0.;
        ds.bc.fisot = 0.;
        ds.bc.albedo = surface_albedo[n];
        c_disort(&ds, &out);
        for (i=0; i<num_levels; ++i)
        {
            uint64_t o = i*(optics->grid.n) + n;
            flux_up[o] = total_solar_irradiance*out.rad[i].flup;
            flux_down[o] = total_solar_irradiance*(out.rad[i].rfldir + out.rad[i].rfldn);
        }
    }

    /* Free allocated memory.*/
    c_disort_out_free(&ds, &out);
    c_disort_state_free(&ds);
    return GRTCODE_SUCCESS;
#else
    char const * mesg = "Program not built with DISORT (returning code %d)."
                        "Please try:\n$ ./configure --enable-disort"
                        " LIBS=<path to libcdisort.a>";
    raise(GRTCODE_COMPILER_ERR, mesg, GRTCODE_COMPILER_ERR);

    /*To quiet compiler warnings.*/
    (void)optics;
    (void)zen_dir;
    (void)surface_albedo;
    (void)total_solar_irradiance;
    (void)solar_flux;
    (void)flux_up;
    (void)flux_down;
#endif
}
