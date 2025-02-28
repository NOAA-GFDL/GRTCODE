#ifndef CLOUDS_LIB_H
#define CLOUDS_LIB_H


int initialize_clouds_lib(char const * beta_path, char const * ice_path,
                          char const * liquid_path, int const * beta_shape);


int finalize_clouds_lib();


int cloud_optics(int const num_bands, double const * band_centers, double const * band_limits,
                 int const num_layers, double const * mean_cloud_fraction,
                 double const * mean_liquid_content, double const *mean_ice_content,
                 double const * overlap, double const liquid_radius,
                 double const * temperature, double * beta_liquid, double * omega_liquid,
                 double * g_liquid, double * beta_ice, double *omega_ice, double * g_ice);


int calculate_overlap(int const num_layers, double const * altitude,
                      double const scale_length, double * alpha);


#endif
