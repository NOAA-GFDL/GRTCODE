#ifndef CLOUDS_LIB_H
#define CLOUDS_LIB_H
#ifndef FP_T_DEFINED
#define FP_T_DEFINED
typedef double fp_t;  /* or use float, as required */
#endif

int initialize_clouds_lib(char const * beta_path, char const * ice_path,
                          char const * liquid_path);


int finalize_clouds_lib();


static double ice_particle_size(double temperature);

int cloud_optics(const double *wavenum,
                 int num_wavenum,
                 int num_layers,
                 const double *mean_cloud_fraction,
                 const double *mean_liquid_content,
                 const double *mean_ice_content,
                 const double *overlap,
                 const double liquid_radius,
                 const double *temperature,
                 double *beta_liquid, double *omega_liquid, double *g_liquid,
                 double *beta_ice, double *omega_ice, double *g_ice);

int calculate_overlap(int const num_layers, double const * altitude,
                      double const scale_length, double * alpha);

#endif
