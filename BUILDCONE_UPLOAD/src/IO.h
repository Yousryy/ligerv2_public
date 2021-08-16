#ifndef _IO_H
#define _IO_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./vardef.h"
float *Read_density_float(char *path,int *Ngrid,float *box_length,float *average_den);
double *Read_density(char *path,int *Ngrid,float *box_length,float *average_den);
double calc_comoving_distance(double red);
struct io_header_1 load_header_gdt(char *fname,int singlefile);
struct survey_functions *Form_surv_struct(char *fname,char delim[]);
struct particle_data *Read_Coordinates(char *fname,int n_gal);
struct magni *Read_magnification(char *path, int n_gal,int r_max);
struct Time_data *Read_dlna(char *path, int n_gal,int r_max);
struct dens_com *Read_nbar_text(char *save_path,int *array_size);
void File_Specs(char *particle_path,int *n_gal,float *r_max);
void Text_file_specs(char *fname,char delim[],int *ncol,int *nrow);
void Save_array(char *save_path,float *array, int array_size,float box_len,float av_den);
void Save_array_text(char *save_path,double *array, int array_size,float grid_size,double *k_bin);
void Read_simulation_parameters(char *music_base);
void Save_nbar_text(char *save_path,struct dens_com *n_bar, int array_size);
void Initiate_spline(struct survey_functions *surv,int nrows);
void Free_spline();
#endif
