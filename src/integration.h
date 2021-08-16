/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn
    Copyright (C) 2021  Mohamed Elkhashab, Argelander Institut für Astronomie, University of Bonn
*/

#ifndef _INTEGRATION_H
#define _INTEGRATION_H
#include "vardef.h"
double Hub(double a,double Om,double OL);
int Find_index_closest_value(struct light_rays_path *path,float chi_particle);
double calc_comoving_distance(double red,double dz,double Om,double OL);
 double simulation_time_interpol(double fieldfarer,double fieldfar,double fieldclose,double fieldcloser,double distfarer,double distfar,double distclose,double distcloser,double dist);
 int get_integration_range(struct light_rays_path *ray,double *rangemin,double *rangemax,double gridsize,int *SnapStart,int *SnapStop);
void integrate_path(struct light_rays_path *ray,double Hub,double *rangemin,double *rangemax,struct potential_container **Pot);
#endif
