#ifndef _CREATE_GALAXY_DENSITY_H
#define _CREATE_GALAXY_DENSITY_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./vardef.h"


void Build_Cone(char *base_liger,char *base_music,char *surv_path,char *save_path,int Ngrid);
float Max_av_shift(int n_gal,struct particle_data *real_par,struct particle_data *shift_par );
float *CIC(struct particle_data *par_struct,struct magni *magn_arr,struct Time_data *T_DLNA,struct survey_functions *srv_funcs,float grid_size,int Ngrid, int n_gal,int mode);
float *Create_galaxy_density(float *dens_real_w,float *dens_shift_w,float *dens_evol,float *dens_mag,struct survey_functions *surv_func,int Ngrid, float gridsize,float max_particle_dist,int mode,float q_value);
double Weight_CIC(float chi_particle,float chi_particle_w,struct magni *mag_arr,struct survey_functions *surv_funcs,int mode);
double get_interpolated_surv_value (int idx_surv,float chi,struct survey_functions *srv_funcs,int n_b_q);
double Cubic_spline(double fieldfarer,double fieldfar,double fieldclose,double fieldcloser,double distfarer,double distfar,double distclose,double distcloser,double dist);
double calc_comoving_distance(double red);
int Find_cell_value(float chi,struct survey_functions *srv_funcs);
int Shift_koor(int idx_x,int idx_y,int idx_z,int Ngrid,int mode);
void Radial_bins(float *galaxy_density,struct survey_functions *srv_funcs,float z_min, float z_max,float grid_size,int Ngrid);
void ADD_shot(float *corr_dens,int Ngrid,float box_length);
#endif
