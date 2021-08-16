//Constants
#ifndef _VARDEF_H
#define _VARDEF_H


#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#define c0 300000
#define PI 3.14159265359
#define hubble_constant (float) 100
#define order 4
#define to_rad     0.01745329252
//Structures
struct particle_data {
  float pos[3];
} ;
struct patricle_gdgt{
  float Pos[3];
  float Vel[3];
  long Id;
};
struct dens_com{
  double n_bar;
  double chi;
  // double *z;
};
struct magni{
  float mag_tot;
  float mag_conv;
  float mag_doppler;
  float mag_tot_obs;
  float mag_doppler_obs;
};
struct Time_data{
  // float time_distance[4];
  float dlna[4];
};
struct survey_functions{
  double redshift;
  double b_lin;
  double chi;
  double n_z;
  double q_mag;
  double e_evl;
} ;
struct survey_functions_splines{
  gsl_interp_accel *acc_red;
  gsl_interp_accel *acc_b_lin;
  gsl_interp_accel *acc_n_z;
  gsl_interp_accel *acc_q_mag;
  gsl_interp_accel *acc_e_bias;
  gsl_spline *spline_red;
  gsl_spline *spline_b_lin;
  gsl_spline *spline_n_z;
  gsl_spline *spline_q_mag;
  gsl_spline *spline_e_bias;
  double *CHI_interp,*Func_CHI;} ;
struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};

//Global Variables

struct Global_VAR {
char liger_path[500];
char sim_path[500];
char surv_func_path[500];
char save_dir_path[500];
char density_directory[200];
float gal_deg, eli_deg;// This was designed for euclid, for other survey please change these numbers
int DM;
float factor_box_size;
float chi_bin_nbar;
int only_dens;
int *dens_to_compute;
int nu_of_cones;
float Q_const;
int Ngrid;
int aliasing;
int SN;
};


struct Global_VAR run;
int over_den_bar;
int Euclid_bins;
double chi_neg,omega_m,omega_L,n_sim_bar;
int surv_number_pts , DM, USESPLINE;
struct survey_functions_splines surv_splines;
#endif
