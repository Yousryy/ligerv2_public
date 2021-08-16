/*  Liger - Calculates the distribution of galaxies (or any targets) on the backward light cone of any observer taking into account leading-order GR effects in the weak-field regime
Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn
Copyright (C) 2021  Mohamed Elkhashab, Argelander Institut für Astronomie, University of Bonn


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

//  This file contains the routines to find the line-of-sight and interpolate the world line of the galaxy with the light cone of the observer

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vardef.h"
#include "potential.h"
#include "integration.h"

int get_grid_cell_to_coordinate(double pos,double gridsize){
	if (pos<0) return -abs(pos/gridsize-0.999999999999); else return abs(pos/gridsize);
}

//Recalculates the distance in case the end position changes due to the coordinate transformation.
void update_ray(struct light_rays_path *path){
	double dx,dy,dz;
	dx = path[0].end[0]-path[0].start[0];
	dy = path[0].end[1]-path[0].start[1];
	dz = path[0].end[2]-path[0].start[2];
	path[0].dist = sqrt(dx*dx+dy*dy+dz*dz);
}

//Sets up the struct light_rays_path for a target defined in "Part"
void initialize_ray(struct light_rays_path *path,struct particle_data Part,double boxsize,double SnapDist,double Hub){
	double dist;

	path[0].TimeDistance = SnapDist;
	path[0].start[0] = run.ObsPos[0];
	path[0].start[1] = run.ObsPos[1];
	path[0].start[2] = run.ObsPos[2];
	path[0].end[0] = PeriDist(Part.Pos[0]-run.ObsPos[0],boxsize);//Particle coordinates are offcentred with respect to the field
	path[0].end[1] = PeriDist(Part.Pos[1]-run.ObsPos[1],boxsize);
	path[0].end[2] = PeriDist(Part.Pos[2]-run.ObsPos[2],boxsize);
	dist=sqrt(path[0].end[0]*path[0].end[0]+path[0].end[1]*path[0].end[1]+path[0].end[2]*path[0].end[2]);

	path[0].normvec[0]=path[0].end[0]/dist;
	path[0].normvec[1]=path[0].end[1]/dist;
	path[0].normvec[2]=path[0].end[2]/dist;
	path[0].vr=(Part.Vel[0]*path[0].normvec[0]+Part.Vel[1]*path[0].normvec[1]+Part.Vel[2]*path[0].normvec[2]);
	path[0].vrHub=path[0].vr/Hub;

	path[0].end[0]+=run.ObsPos[0];
	path[0].end[1]+=run.ObsPos[1];//Put to proper center with respect to grav potential
	path[0].end[2]+=run.ObsPos[2];
	path[0].dist=dist;
	path[0].dlna =0;
	}

//Finds all cells, and the lenght of the intersecting line for each cell, between the start and end position defined in the stuct light_rays_path
void find_cells_along_ray(struct light_rays_path *path,double gridsize,int Ngrid){
	int gs[3],ge[3],c,gx,gy,gz;
	double stepx,stepy,stepz,lenx,leny,lenz,norm[3];
	gs[0]=get_grid_cell_to_coordinate(path[0].start[0],gridsize);
	gs[1]=get_grid_cell_to_coordinate(path[0].start[1],gridsize);
	gs[2]=get_grid_cell_to_coordinate(path[0].start[2],gridsize);
	ge[0]=get_grid_cell_to_coordinate(path[0].end[0],gridsize);
	ge[1]=get_grid_cell_to_coordinate(path[0].end[1],gridsize);
	ge[2]=get_grid_cell_to_coordinate(path[0].end[2],gridsize);

	norm[0]=path[0].end[0]-path[0].start[0];
	norm[1]=path[0].end[1]-path[0].start[1];
	norm[2]=path[0].end[2]-path[0].start[2];
	lenx=sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
	norm[0]/=lenx;
	norm[1]/=lenx;
	norm[2]/=lenx;
	stepx=gridsize*sqrt(1.0+norm[1]*norm[1]/(norm[0]*norm[0])+norm[2]*norm[2]/(norm[0]*norm[0]));
	stepy=gridsize*sqrt(1.0+norm[0]*norm[0]/(norm[1]*norm[1])+norm[2]*norm[2]/(norm[1]*norm[1]));
	stepz=gridsize*sqrt(1.0+norm[0]*norm[0]/(norm[2]*norm[2])+norm[1]*norm[1]/(norm[2]*norm[2]));

	if (finite(stepx)==1) gx=norm[0]/fabs(norm[0]); else gx=0;
	if (finite(stepy)==1) gy=norm[1]/fabs(norm[1]); else gy=0;
	if (finite(stepz)==1) gz=norm[2]/fabs(norm[2]); else gz=0;

	if (PeriDist(path[0].end[0],gridsize*Ngrid)/gridsize==ge[0]){
		if (gx<0) ge[0]--;
		else if ((gx>0)&&(ge[0]<gs[0])) ge[0]++;
	}
	if (PeriDist(path[0].end[1],gridsize*Ngrid)/gridsize==ge[1]){
		if (gy<0) ge[1]--;
		else if ((gy>0)&&(ge[1]<gs[1])) ge[1]++;
	}
	if (PeriDist(path[0].end[2],gridsize*Ngrid)/gridsize==ge[2]){
		if (gz<0) ge[2]--;
		else if ((gz>0)&&(ge[2]<gs[2])) ge[2]++;
	}
	path[0].CellsAlongPath=abs(gs[0]-ge[0])+abs(gs[1]-ge[1])+abs(gs[2]-ge[2])+1; //+1 for the final cell

	if (!(path[0].PathCell=malloc(path[0].CellsAlongPath*sizeof(int)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
	if (!(path[0].PathIntersec=malloc(path[0].CellsAlongPath*sizeof(float)))) {printf("\nAllocation of memory failed.\n"); exit(0);}

	lenx=fabs(path[0].start[0]-(gs[0]+(gx+1.)/2)*gridsize)*stepx/gridsize;
	leny=fabs(path[0].start[1]-(gs[1]+(gy+1.)/2)*gridsize)*stepy/gridsize;
	lenz=fabs(path[0].start[2]-(gs[2]+(gz+1.)/2)*gridsize)*stepz/gridsize;
	if (finite(lenx)!=1.) lenx=path[0].dist*10;
	if (finite(leny)!=1.) leny=path[0].dist*10;//Something which is finite but large enough.
	if (finite(lenz)!=1.) lenz=path[0].dist*10;

	c=0;
	while ((gs[0]!=ge[0])||(gs[1]!=ge[1])||(gs[2]!=ge[2])) {
		if (c>=path[0].CellsAlongPath) {printf("Failed to reach end in the ray tracing!\nStart:	%lf	%lf	%lf\nEnd:	%lf	%lf	%lf\nEst. #cells: %i\nSteplenght:	%lf	%lf	%lf\nGridCurr:	%i	%i	%i\nGridDirect:	%i	%i	%i\nFinal Grid Pos:	%lf	%lf	%lf\nExact?:	%i	%i	%i\n",path[0].start[0],path[0].start[1],path[0].start[2],path[0].end[0],path[0].end[1],path[0].end[2],path[0].CellsAlongPath,stepx,stepy,stepz,gs[0],gs[1],gs[2],gx,gy,gz,path[0].end[0]/gridsize,path[0].end[1]/gridsize,path[0].end[2]/gridsize,ge[0]==path[0].end[0]/gridsize,ge[1]==path[0].end[1]/gridsize,ge[2]==path[0].end[2]/gridsize); fflush(stdout); exit(0);}
		if (lenx < leny) {
			if (lenx < lenz) {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=lenx;
				gs[0]+=gx;
				lenx+=stepx;
				c++;
			} else if (lenx > lenz) {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=lenz;
				gs[2]+=gz;
				lenz+=stepz;
				c++;
			} else {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=lenz;
				gs[0]+=gx;
				gs[2]+=gz;
				lenx+=stepx;
				lenz+=stepz;
				c++;
			}
		} else if ((lenx != leny)||(leny > lenz)) {
			if (leny < lenz) {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=leny;
				gs[1]+=gy;
				leny+=stepy;
				c++;
			} else if (leny > lenz) {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=lenz;
				gs[2]+=gz;
				lenz+=stepz;
				c++;
			} else {
				path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
				path[0].PathIntersec[c]=lenz;
				gs[1]+=gy;
				gs[2]+=gz;
				leny+=stepy;
				lenz+=stepz;
				c++;
			}
		} else if (lenx == lenz) {
			path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
			path[0].PathIntersec[c]=lenz;
			gs[0]+=gx;
			gs[1]+=gy;
			gs[2]+=gz;
			lenx+=stepx;
			leny+=stepy;
			lenz+=stepz;
			c++;
		} else {
			path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
			path[0].PathIntersec[c]=lenx;
			gs[0]+=gx;
			gs[1]+=gy;
			lenx+=stepx;
			leny+=stepy;
			c++;
		}
	}
	if (path[0].PathIntersec[c-1]<path[0].dist){
		if (c>=path[0].CellsAlongPath) {printf("Overshooting path array!\nStart:	%lf	%lf	%lf\nEnd:	%lf	%lf	%lf\nEst. #cells: %i\nSteplenght:	%lf	%lf	%lf\nGridCurr:	%i	%i	%i\nGridDirect:	%i	%i	%i\nFinal Grid Pos:	%lf	%lf	%lf\nExact?:	%i	%i	%i\n",path[0].start[0],path[0].start[1],path[0].start[2],path[0].end[0],path[0].end[1],path[0].end[2],path[0].CellsAlongPath,stepx,stepy,stepz,gs[0],gs[1],gs[2],gx,gy,gz,path[0].end[0]/gridsize,path[0].end[1]/gridsize,path[0].end[2]/gridsize,ge[0]==path[0].end[0]/gridsize,ge[1]==path[0].end[1]/gridsize,ge[2]==path[0].end[2]/gridsize); fflush(stdout); exit(0);}
		path[0].PathCell[c]=shift_koor(gs[0],gs[1],gs[2],Ngrid);
		path[0].PathIntersec[c]=path[0].dist;
		c++;
	}
	path[0].CellsAlongPath=c;
	path[0].PathIntersec_max_line = path[0].PathIntersec[path[0].CellsAlongPath-1];
}

//Adds or removes (bit of a misnomer honestly) observor velocity shifts post integration
//note that timedistance is already in distance units
void add_obs_vel(double hub_value,struct light_rays_path *ray,int ifRel,int ifadd) {
	double MaxLineOfSight	= ray[0].PathIntersec_max_line;
	double vr_obs = (run.ObsVel[0]*ray[0].normvec[0]+run.ObsVel[1]*ray[0].normvec[1]+run.ObsVel[2]*ray[0].normvec[2]);
	double potential_norm = run.Pot_Obs*SPEEDLIGHTSQUAREINV;
	if (ifRel == 1){
		if (ifadd == 1){
			ray[0].end[0] -= ray[0].normvec[0]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[0] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[1] -= ray[0].normvec[1]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[1] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[2] -= ray[0].normvec[2]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[2] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].TimeDistance -= SPEEDLIGHT*potential_norm/hub_value - vr_obs/hub_value;
			ray[0].dlna +=  potential_norm - vr_obs/SPEEDLIGHT;
		}else{
			ray[0].end[0] += ray[0].normvec[0]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[0] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[1] += ray[0].normvec[1]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[1] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[2] += ray[0].normvec[2]*((MaxLineOfSight/SPEEDLIGHT + 1./hub_value)*vr_obs - (2.*MaxLineOfSight + SPEEDLIGHT/hub_value)*potential_norm) - run.ObsVel[2] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].TimeDistance += SPEEDLIGHT*potential_norm/hub_value - vr_obs/hub_value;
			ray[0].dlna -=  potential_norm - vr_obs/SPEEDLIGHT;
		}
	}else{
		if (ifadd == 1){
			ray[0].end[0] -= ray[0].normvec[0]*vr_obs/hub_value - run.ObsVel[0] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[1] -= ray[0].normvec[1]*vr_obs/hub_value - run.ObsVel[1] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[2] -= ray[0].normvec[2]*vr_obs/hub_value - run.ObsVel[2] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].TimeDistance -= - vr_obs/(hub_value);
			ray[0].dlna +=  - vr_obs/SPEEDLIGHT;
		}
		else{
			ray[0].end[0] += ray[0].normvec[0]*vr_obs/hub_value - run.ObsVel[0] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[1] += ray[0].normvec[1]*vr_obs/hub_value - run.ObsVel[1] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].end[2] += ray[0].normvec[2]*vr_obs/hub_value - run.ObsVel[2] * MaxLineOfSight/SPEEDLIGHT;
			ray[0].TimeDistance += - vr_obs/(hub_value);
			ray[0].dlna -=  - vr_obs/SPEEDLIGHT;
		}
	}
}
//Finds the position (and magnification) of the galaxy defined in ray[0] on the light cone by interpolating ds=dr^2-dt^2 to ds=0
//Mode, 0:Relativistic with vel_obs,1:Relativistic only,2:Doppler shift with observer,3:Doppler shift only
int light_cone_intersection(struct light_rays_path *ray,struct header_snapshot Head,struct light_cone_object *Result,int Mode){
	int b,SnapLC,SnapLCm2,SnapLCm1,SnapLCp1;
	double dist,Mindist,InterpolWL[4];

	SnapLC=-10;
	Mindist=-Head.BoxSize*Head.BoxSize*9;
	for (b=0;b<run.NSnap;++b) {// Find between which snapshots the world line and the light cone intersect
		dist=ray[b].dist*ray[b].dist-ray[b].TimeDistance*ray[b].TimeDistance;
		if ((dist>Mindist)&&(dist<0.0)) {Mindist=dist; SnapLC=b;}
	}

	if (SnapLC<0) return -1; //World line too far away, no time-like snapshot.
	// if ((SnapLC==0)&&(ray[SnapLC].TimeDistance*ray[SnapLC].TimeDistance > ray[SnapLC].dist*ray[SnapLC].dist)){
	// 	printf("\b\b\b\bWarning: particle to close to observer that shift is bigger than dist\n000");
	// 	return -1; //particle to close to observer that shift is bigger than dist.
	// }
	if ((SnapLC==0)&&(ray[0].TimeDistance*ray[0].TimeDistance > ray[0].dist*ray[0].dist)){
		//printf("\b\b\b\bWarning: Particle track ended before crossing light cone!\n000");
		//Ended_before_crossing_light_cone++;
		Result[0].SimPos[0]=PeriDist(ray[0].end[0]-ray[0].start[0],Head.BoxSize);
		Result[0].SimPos[1]=PeriDist(ray[0].end[1]-ray[0].start[1],Head.BoxSize);
		Result[0].SimPos[2]=PeriDist(ray[0].end[2]-ray[0].start[2],Head.BoxSize);
		Result[0].Magni[0]=ray[0].Magni[0];
		Result[0].Magni[1]=ray[0].Magni[1];
		Result[0].Magni[2]=ray[0].Magni[2];
		return 2;
	}
	if ((SnapLC==0)||(SnapLC>=run.NSnap)) {
		printf("Error: Out of range snapshot supscript! Rel %i\n",SnapLC);
		for (b=0;b<run.NSnap;b++) printf("%.1e ",ray[b].dist*ray[b].dist-ray[b].TimeDistance*ray[b].TimeDistance);
		printf("\n");
		fflush(stdout);
		exit(0);
	}
	// SnapLC,SnapLCm2,SnapLCm1,SnapLCp1;
	SnapLCm1 = SnapLC-1;
	if (SnapLC==1) SnapLCm2=SnapLC-1; else SnapLCm2=SnapLC-2;
	if (SnapLC==run.NSnap-1) SnapLCp1=SnapLC; else SnapLCp1=SnapLC+1;

	InterpolWL[0]=ray[SnapLCm2].dist*ray[SnapLCm2].dist-ray[SnapLCm2].TimeDistance*ray[SnapLCm2].TimeDistance;
	InterpolWL[1]=ray[SnapLCm1].dist*ray[SnapLCm1].dist-ray[SnapLCm1].TimeDistance*ray[SnapLCm1].TimeDistance;
	InterpolWL[2]=ray[SnapLC].dist*ray[SnapLC].dist-ray[SnapLC].TimeDistance*ray[SnapLC].TimeDistance;
	InterpolWL[3]=ray[SnapLCp1].dist*ray[SnapLCp1].dist-ray[SnapLCp1].TimeDistance*ray[SnapLCp1].TimeDistance;

	Result[0].SimPos[0]=PeriDist(simulation_time_interpol(ray[SnapLCp1].end[0],ray[SnapLC].end[0],ray[SnapLCm1].end[0],ray[SnapLCm2].end[0],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0)-run.ObsPos[0],Head.BoxSize);
	Result[0].SimPos[1]=PeriDist(simulation_time_interpol(ray[SnapLCp1].end[1],ray[SnapLC].end[1],ray[SnapLCm1].end[1],ray[SnapLCm2].end[1],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0)-run.ObsPos[1],Head.BoxSize);
	Result[0].SimPos[2]=PeriDist(simulation_time_interpol(ray[SnapLCp1].end[2],ray[SnapLC].end[2],ray[SnapLCm1].end[2],ray[SnapLCm2].end[2],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0)-run.ObsPos[2],Head.BoxSize);


	//Interpolation of timedistance and dln a
	Result[0].TimePosition[0]=simulation_time_interpol(ray[SnapLCp1].TimeDistance,ray[SnapLC].TimeDistance,ray[SnapLCm1].TimeDistance,ray[SnapLCm2].TimeDistance,InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
	Result[0].TimePosition[1]=simulation_time_interpol(ray[SnapLCp1].dlna,ray[SnapLC].dlna,ray[SnapLCm1].dlna,ray[SnapLCm2].dlna,InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);


	// 	if ((ray[SnapLCm2].D3shift[0]==0.0)||(ray[SnapLCm1].D3shift[0]==0.0)||(ray[SnapLC].D3shift[0]==0.0)||(ray[SnapLCp1].D3shift[0]==0.0)) {		printf("Warning No shifts computed for particle!\n");
	// 	printf("%i	%le	%le	%le	%le	%le %le %le %le\n",SnapLCm2,ray[SnapLCm2].DTshift,ray[SnapLCm2].DRshift,ray[SnapLCm2].D3shift[0],ray[SnapLCm2].D3shift[1],ray[SnapLCm2].D3shift[2],ray[SnapLCm2].dist,ray[SnapLCm2].TimeDistance,ray[SnapLCm2].end[2]);
	// 	printf("%i	%le	%le	%le	%le	%le %le %le %le\n",SnapLCm1,ray[SnapLCm1].DTshift,ray[SnapLCm1].DRshift,ray[SnapLCm1].D3shift[0],ray[SnapLCm1].D3shift[1],ray[SnapLCm1].D3shift[2],ray[SnapLCm1].dist,ray[SnapLCm1].TimeDistance,ray[SnapLCm1].end[2]);
	// 	printf("%i	%le	%le	%le	%le	%le %le %le %le\n",SnapLC,ray[SnapLC].DTshift,ray[SnapLC].DRshift,ray[SnapLC].D3shift[0],ray[SnapLC].D3shift[1],ray[SnapLC].D3shift[2],ray[SnapLC].dist,ray[SnapLC].TimeDistance,ray[SnapLC].end[2]);
	// 	printf("%i	%le	%le	%le	%le	%le %le %le %le\n",SnapLCp1,ray[SnapLCp1].DTshift,ray[SnapLCp1].DRshift,ray[SnapLCp1].D3shift[0],ray[SnapLCp1].D3shift[1],ray[SnapLCp1].D3shift[2],ray[SnapLCp1].dist,ray[SnapLCp1].TimeDistance,ray[SnapLCp1].end[2]);
	// 	fflush(stdout);
	// 	if ((ray[SnapLCm2].dist>run.mindist)||(ray[SnapLCm1].dist>run.mindist)||(ray[SnapLC].dist>run.mindist)||(ray[SnapLCp1].dist>run.mindist)) exit(0);
	// 	else {printf("Outside range of cone %f\t%f\n",run.mindist,run.maxdist);fflush(stdout);}
	// }
	//Internal safety conditions if shifts are not computed the data of the particle will be provided in the stdout but the code won't crash
	//replace if commands with previous code if you want code to shut down for any particle with no integral
	//Mode, 0:Relativistic with vel_obs,1:Relativistic only,2:Doppler shift with observer,3:Doppler shift only
	if (Mode==1) {
		if ((ray[SnapLCm2].D3shift[0]==0.0)||(ray[SnapLCm1].D3shift[0]==0.0)||(ray[SnapLC].D3shift[0]==0.0)||(ray[SnapLCp1].D3shift[0]==0.0)) {
			return 100;
		}
		Result[0].Magni[0]=simulation_time_interpol(ray[SnapLCp1].Magni[0],ray[SnapLC].Magni[0],ray[SnapLCm1].Magni[0],ray[SnapLCm2].Magni[0],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
		Result[0].Magni[1]=simulation_time_interpol(ray[SnapLCp1].Magni[1],ray[SnapLC].Magni[1],ray[SnapLCm1].Magni[1],ray[SnapLCm2].Magni[1],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
		Result[0].Magni[2]=simulation_time_interpol(ray[SnapLCp1].Magni[2],ray[SnapLC].Magni[2],ray[SnapLCm1].Magni[2],ray[SnapLCm2].Magni[2],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
	}else if (Mode==0) {
		if ((ray[SnapLCm2].D3shift[0]==0.0)||(ray[SnapLCm1].D3shift[0]==0.0)||(ray[SnapLC].D3shift[0]==0.0)||(ray[SnapLCp1].D3shift[0]==0.0)) {
			return 100;
		}
		Result[0].Magni[0]=simulation_time_interpol(ray[SnapLCp1].Magni[3],ray[SnapLC].Magni[3],ray[SnapLCm1].Magni[3],ray[SnapLCm2].Magni[3],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
		Result[0].Magni[1]=simulation_time_interpol(ray[SnapLCp1].Magni[4],ray[SnapLC].Magni[4],ray[SnapLCm1].Magni[4],ray[SnapLCm2].Magni[4],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
		Result[0].Magni[2]=simulation_time_interpol(ray[SnapLCp1].Magni[5],ray[SnapLC].Magni[5],ray[SnapLCm1].Magni[5],ray[SnapLCm2].Magni[5],InterpolWL[3],InterpolWL[2],InterpolWL[1],InterpolWL[0],0.0);
	}
	dist=sqrt(Result[0].SimPos[0]*Result[0].SimPos[0]+Result[0].SimPos[1]*Result[0].SimPos[1]+Result[0].SimPos[2]*Result[0].SimPos[2]);
  if ((dist<run.maxdist)&&(dist>run.mindist)) return 0;
	else return 1;
}

//A wrapper around "light_cone_intersection" to find the real- and redshift-space light-cone intersection.
int get_all_light_cone_intersection(struct light_rays_path *ray,struct header_snapshot Head,struct light_cone_object *Result){
	int b,ReSim,ReVel,ReRel,ReRel_obs,ReVel_obs;
	struct light_cone_object Step;
	//Mode = Relativistic with observor results
	ReRel_obs = light_cone_intersection(ray,Head,&Step,0);
	if (ReRel_obs>=0) {
		Result[0].RelPos_obs[0]=Step.SimPos[0];
		Result[0].RelPos_obs[1]=Step.SimPos[1];
		Result[0].RelPos_obs[2]=Step.SimPos[2];
		Result[0].Magni[3]=Step.Magni[0];
		Result[0].Magni[4]=Step.Magni[1];
		Result[0].Magni[5]=Step.Magni[2];
		Result[0].TimePosition[3] = Step.TimePosition[0];
		Result[0].TimePosition[7] = Step.TimePosition[1];
	}else return -1;
	for (b=0;b<run.NSnap;++b){
		if(ray[b].D3shift[0]==0.0) continue;
		add_obs_vel(Hubarray[b],&ray[b],1,0);
		update_ray(&ray[b]);
	}
	//Mode = Relativistic without observer results

	ReRel = light_cone_intersection(ray,Head,&Step,1);
	if (ReRel>=0) {
		Result[0].RelPos[0]=Step.SimPos[0];
		Result[0].RelPos[1]=Step.SimPos[1];
		Result[0].RelPos[2]=Step.SimPos[2];
		Result[0].Magni[0]=Step.Magni[0];
		Result[0].Magni[1]=Step.Magni[1];
		Result[0].Magni[2]=Step.Magni[2];
		Result[0].TimePosition[2] = Step.TimePosition[0];
		Result[0].TimePosition[6] = Step.TimePosition[1];
	}else return -1;
	//remove GR shifts and add only observer velocity effects
	for (b=0;b<run.NSnap;++b){
		if(ray[b].D3shift[0]==0.0) continue;
		ray[b].end[0]+=ray[b].DRshift*ray[b].normvec[0]+ray[b].D3shift[0]; //Shift position back.
		ray[b].end[1]+=ray[b].DRshift*ray[b].normvec[1]+ray[b].D3shift[1];
		ray[b].end[2]+=ray[b].DRshift*ray[b].normvec[2]+ray[b].D3shift[2];
		ray[b].TimeDistance+=ray[b].DTshift;
		ray[b].dlna-=ray[b].DTshift*Hubarray[b]/SPEEDLIGHT;
		add_obs_vel(Hubarray[b],&ray[b],0,1);
		update_ray(&ray[b]);
	}

	//Mode = velocity with observer velocity
	ReVel_obs=light_cone_intersection(ray,Head,&Step,2);
	if (ReVel_obs>=0) {
		Result[0].VelPos_obs[0]=Step.SimPos[0];
		Result[0].VelPos_obs[1]=Step.SimPos[1];
		Result[0].VelPos_obs[2]=Step.SimPos[2];
		Result[0].TimePosition[1] = Step.TimePosition[0];
		Result[0].TimePosition[5] = Step.TimePosition[1];
	}else return -1;
	//remove observer velocity
	for (b=0;b<run.NSnap;++b){
		if(ray[b].D3shift[0]==0.0) continue;
		add_obs_vel(Hubarray[b],&ray[b],0,0);
		update_ray(&ray[b]);

	}
	//Mode =Velocity results
	ReVel=light_cone_intersection(ray,Head,&Step,3);
	if (ReVel>=0) {
		Result[0].VelPos[0]=Step.SimPos[0];
		Result[0].VelPos[1]=Step.SimPos[1];
		Result[0].VelPos[2]=Step.SimPos[2];
		Result[0].TimePosition[0] = Step.TimePosition[0];
		Result[0].TimePosition[4] = Step.TimePosition[1];
	}else return -1;
//remove velocity effects
	for (b=0;b<run.NSnap;++b){
		ray[b].end[0]-=ray[b].vrHub*ray[b].normvec[0]; //Shift position back.
		ray[b].end[1]-=ray[b].vrHub*ray[b].normvec[1];
		ray[b].end[2]-=ray[b].vrHub*ray[b].normvec[2];
		ray[b].TimeDistance+=ray[b].vrHub;
		ray[b].dlna -= ray[b].vrHub *Hubarray[b]/SPEEDLIGHT;
		if ((fabs(ray[b].dlna)>1e-3)&&(fabs(ray[b].D3shift[0])>1e-3)) {printf("WARNING: DLNA NOT PROPERLY REMOVED%e \n",ray[b].dlna);}
		update_ray(&ray[b]);
	}
	// Mode = simulation results (configuration space)
	ReSim=light_cone_intersection(ray,Head,&Step,4);


	if (ReSim>=0) {
		Result[0].SimPos[0]=Step.SimPos[0];
		Result[0].SimPos[1]=Step.SimPos[1];
		Result[0].SimPos[2]=Step.SimPos[2];
	} else return -1;
	//add condition particles too close
	if ((ReRel==100)||(ReVel==100)||(ReSim==100)||(ReRel_obs==100)||(ReVel_obs==100)) return 100;
	else if ((ReRel==0)||(ReVel==0)||(ReSim==0)||(ReRel_obs==0)||(ReVel_obs==0)) return 0;
	else return 1;
}
