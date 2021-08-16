/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn
    Copyright (C) 2021  Mohamed Elkhashab, Argelander Institut für Astronomie, University of Bonn
*/

#ifndef _LIGHT_RAYS_H
#define _LIGHT_RAYS_H
void update_ray(struct light_rays_path *path);
void initialize_ray(struct light_rays_path *path,struct particle_data Part,double boxsize,double SnapDist,double Hub);
void find_cells_along_ray(struct light_rays_path *path,double gridsize,int Ngrid);
int light_cone_intersection(struct light_rays_path *ray,struct header_snapshot Head,struct light_cone_object *Result,int ifMagni);
int get_all_light_cone_intersection(struct light_rays_path *ray,struct header_snapshot Head,struct light_cone_object *Result);
void add_obs_vel(double hub_value,struct light_rays_path *ray,int ifRel,int ifadd);
#endif
