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

//  This file contains the routines to perform the integration along the line of sight.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./vardef.h"

double Hub(double a,double Om,double OL){
	/*Calculates Hubble(a) but he puts H_0 to 100.0 bby construction*/
return 100.0*sqrt(Om/a/a/a+(1.0-Om-OL)/a/a+OL);
//hub=(0.0001022729969);// in units of 1/Myr
}
int Find_index_closest_value(struct light_rays_path *path,float chi_particle){
  float chi,chi_first,chi_last,chi_middle;
  int first=0,last,middle,idx_of_chi ,old_mid;
  chi =  chi_particle;
  last=path[0].CellsAlongPath;
  if (path[0].PathIntersec[last]<chi_particle) return last-1;
  middle = (first+last)/2;
  chi_middle =path[0].PathIntersec[middle];
  while (old_mid!=middle) {
    old_mid = middle;
    if (chi>chi_middle) {
      first = middle;
      middle = (first+last)/2;
      chi_middle =path[0].PathIntersec[middle];
    }else if (chi<chi_middle){
      last = middle;
      middle = (first+last)/2;
      chi_middle =path[0].PathIntersec[middle];
    }
  }
  idx_of_chi= old_mid;
  while(path[0].PathIntersec[idx_of_chi]>chi) idx_of_chi--;
  return idx_of_chi;
}
double calc_comoving_distance(double red,double dz,double Om,double OL){
int i=0;
double z,result=0.0;
result=pow(Om*pow(red+1.,3)+(red+1.)*(red+1.)*(1.-Om-OL)+OL,-0.5)/2.0+0.5;
for (z=dz;z<red-dz*0.1;z+=dz,i++) result+=pow(Om*pow(z+1.,3)+(z+1.)*(z+1.)*(1.-Om-OL)+OL,-0.5);
return result*dz*2997.92458; // Mpc/h
}

//interpolates four points to the desired location.
// switched input values
 // double simulation_time_interpol(double fieldfarer,double fieldfar,double fieldclose,double fieldcloser,double distfarer,double distfar,double distclose,double distcloser,double dist){
  double simulation_time_interpol(double fieldcloser,double fieldclose,double fieldfar,double fieldfarer,double distcloser,double distclose,double distfar,double distfarer,double dist){
	double pg,g,p,c,d,dfdxc,dfdxf;
// g=distfarer-distcloser;
g=distfar-distclose;
p=(dist-distclose);
if ((p < 0)||(p > g)) {
	#pragma omp critical
	{
	printf("Error, out of range interpolation!\nWanted:	%le	(%le	%le)\nBoundary:	%le	%le	%le	%le\nField:	%le	%le	%le	%le\n",dist,p,g,distcloser,distclose,distfar,distfarer,fieldcloser,fieldclose,fieldfar,fieldfarer);
	fflush(stdout);
	exit(0);
	}
}
if (fieldfarer == fieldfar) dfdxf=(fieldfar-fieldclose)/g;
else {
	c=g/(distfarer-distfar);
	dfdxf=(fieldclose-fieldfar-(fieldfarer-fieldfar)*c*c)/(-g*(1.0+c));
}
if (fieldcloser == fieldclose) dfdxc=(fieldfar-fieldclose)/g;
else {
	c=(distcloser-distclose)/g;
	dfdxc=(fieldcloser-fieldclose-(fieldfar-fieldclose)*c*c)/((distcloser-distclose)*(1.0-c));
}
pg=p/g;
c=(3*(fieldfar-fieldclose)-dfdxf*g-2.*dfdxc*g);
d=(dfdxf*g+dfdxc*g-2*(fieldfar-fieldclose));
return fieldclose+dfdxc*p+c*pg*pg+d*pg*pg*pg;
}

//To minimize the numer of integrations, which have to be performed, this function estimates the four snapshots around the crossing of the world line and the observer's light cone.
int get_integration_range(struct light_rays_path *ray,double *rangemin,double *rangemax,double gridsize,int *SnapStart,int *SnapStop){
int i=-1,b,Start,Stop;
double dist,Mindist=-pow(gridsize*run.Ngrid,2.)*9.;
for (b=0;b<run.NSnap;++b) {// Find between which snapshots the world line and the light cone intersect
	dist=ray[b].dist*ray[b].dist-ray[b].TimeDistance*ray[b].TimeDistance;
	if ((dist>Mindist)&&(dist<0.0)) {Mindist=dist; i=b;}
}
if (i==-1) return 1;
if ((ray[i].dist-fabs(ray[i].vrHub) > run.maxdist) || (ray[i].dist+fabs(ray[i].vrHub) < run.mindist)) return 2;
Start=i-2;
Stop=i+3;
if (i>0) if (pow(ray[i].dist+ray[i].vrHub,2.)-pow(ray[i].TimeDistance-ray[i].vrHub,2.) < 0) --Start;
if (i<run.NSnap -1) if (pow(ray[i+1].dist+ray[i+1].vrHub,2.)-pow(ray[i+1].TimeDistance-ray[i+1].vrHub,2.) > 0) ++Stop;
if (Start < 0) Start=0;
if (Stop > run.NSnap) Stop=run.NSnap;
SnapStart[0]=Start;
SnapStop[0]=Stop;
return 0;
}

//This function integrates the path, which was defined by "find_cells_along_ray" in light_rays.c

void integrate_path(struct light_rays_path *ray,double Hub,double *rangemin,double *rangemax,struct potential_container **Pot){
double MaxLineOfSight,IntegratedPot[12],IntStep[10],dr,Rcent,sum ,vr_obs,vH_obs;
int idx_0=0,idx_1=0,idx_2=1,idx_3=2; // smaller number indicates closer snapshot
int n,l,start_loop=0;
MaxLineOfSight=ray[0].PathIntersec[ray[0].CellsAlongPath-1];
ray[0].DRshift = (ray[0].pot-run.Pot_Obs)/(SPEEDLIGHT*Hub) - 2*MaxLineOfSight *run.Pot_Obs* SPEEDLIGHTSQUAREINV ;
ray[0].DTshift = -1.*(ray[0].pot-run.Pot_Obs)/(SPEEDLIGHT*Hub);//added \phi_obs
ray[0].Magni[0]=1.+2.*ray[0].pot*SPEEDLIGHTSQUAREINV;//add \phi_0
vr_obs = (run.ObsVel[0]*ray[0].normvec[0]+run.ObsVel[1]*ray[0].normvec[1]+run.ObsVel[2]*ray[0].normvec[2]);
for (int i=0;i<11;i++) IntegratedPot[i]=0.0;
idx_0=0;
if (rangemax[0]>0.0){while(rangemax[0]>ray[0].PathIntersec[start_loop]) start_loop++;}// In case first snashot z!=0.0
for (n=start_loop;(n<ray[0].CellsAlongPath)&&(idx_2<run.NSnap);n++){
  l=ray[0].PathCell[n];
  if (n==0) dr=ray[0].PathIntersec[n]; else dr=ray[0].PathIntersec[n]-ray[0].PathIntersec[n-1];
  if (n==0) Rcent=ray[0].PathIntersec[n]; else Rcent=(ray[0].PathIntersec[n]+ray[0].PathIntersec[n-1])*0.5;
  if(rangemax[idx_2] < Rcent){
    idx_1+=1;
    idx_2+=1;
  }
  if (idx_1==0) idx_0=idx_1; else idx_0=idx_1-1;
  if (idx_2==run.NSnap-1) idx_3=idx_2; else idx_3=idx_2+1;
  if (idx_2==run.NSnap) break;
  if((Rcent>rangemax[idx_2])||(Rcent<rangemax[idx_1]))
  {printf("ERROR IN INTEGRAL LOOP,%f,%f,%f",Rcent,rangemax[idx_1],rangemax[idx_2]);
  exit(0);}
  IntegratedPot[0]+=simulation_time_interpol(Pot[l][idx_0].val,Pot[l][idx_1].val,Pot[l][idx_2].val,Pot[l][idx_3].val,rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[0]=simulation_time_interpol(Pot[l][idx_0].timeD,Pot[l][idx_1].timeD,Pot[l][idx_2].timeD,Pot[l][idx_3].timeD,rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntegratedPot[1]+=IntStep[0];
  IntegratedPot[2]+=IntStep[0]*Rcent;
  IntStep[1]=simulation_time_interpol(Pot[l][idx_0].firD[0],Pot[l][idx_1].firD[0],Pot[l][idx_2].firD[0],Pot[l][idx_3].firD[0],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntegratedPot[3]+=IntStep[1];
  IntegratedPot[4]+=IntStep[1]*Rcent;
  IntStep[2]=simulation_time_interpol(Pot[l][idx_0].firD[1],Pot[l][idx_1].firD[1],Pot[l][idx_2].firD[1],Pot[l][idx_3].firD[1],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntegratedPot[5]+=IntStep[2];
  IntegratedPot[6]+=IntStep[2]*Rcent;

  IntStep[3]=simulation_time_interpol(Pot[l][idx_0].firD[2],Pot[l][idx_1].firD[2],Pot[l][idx_2].firD[2],Pot[l][idx_3].firD[2],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntegratedPot[7]+=IntStep[3];
  IntegratedPot[8]+=IntStep[3]*Rcent;

  IntStep[4]=simulation_time_interpol(Pot[l][idx_0].secD[0],Pot[l][idx_1].secD[0],Pot[l][idx_2].secD[0],Pot[l][idx_3].secD[0],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[5]=simulation_time_interpol(Pot[l][idx_0].secD[1],Pot[l][idx_1].secD[1],Pot[l][idx_2].secD[1],Pot[l][idx_3].secD[1],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[6]=simulation_time_interpol(Pot[l][idx_0].secD[2],Pot[l][idx_1].secD[2],Pot[l][idx_2].secD[2],Pot[l][idx_3].secD[2],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[7]=simulation_time_interpol(Pot[l][idx_0].croD[0],Pot[l][idx_1].croD[0],Pot[l][idx_2].croD[0],Pot[l][idx_3].croD[0],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[8]=simulation_time_interpol(Pot[l][idx_0].croD[1],Pot[l][idx_1].croD[1],Pot[l][idx_2].croD[1],Pot[l][idx_3].croD[1],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  IntStep[9]=simulation_time_interpol(Pot[l][idx_0].croD[2],Pot[l][idx_1].croD[2],Pot[l][idx_2].croD[2],Pot[l][idx_3].croD[2],rangemax[idx_0],rangemax[idx_1],rangemax[idx_2],rangemax[idx_3],Rcent)*dr;
  sum=IntStep[4]+IntStep[5]+IntStep[6];
  sum+=-2./MaxLineOfSight*(ray[0].normvec[0]*IntStep[1]+ray[0].normvec[1]*IntStep[2]+ray[0].normvec[2]*IntStep[3]);
  sum-=ray[0].normvec[0]*ray[0].normvec[0]*IntStep[4] + ray[0].normvec[1]*ray[0].normvec[1]*IntStep[5] + ray[0].normvec[2]*ray[0].normvec[2]*IntStep[6] + 2.*ray[0].normvec[0]*ray[0].normvec[1]*IntStep[7] + 2.*ray[0].normvec[0]*ray[0].normvec[2]*IntStep[8] + 2.*ray[0].normvec[1]*ray[0].normvec[2]*IntStep[9];
  IntegratedPot[9]+=sum*Rcent;
  IntegratedPot[10]+=sum*Rcent*Rcent;
}

//Now we compute the coordinate shift from the integrated quantities.
ray[0].DRshift+=4.0*IntegratedPot[0]*SPEEDLIGHTSQUAREINV;//why 4 because when you add it to D3 shift you get four
ray[0].DRshift+=2.0*(1.0/(SPEEDLIGHT*Hub)+MaxLineOfSight*SPEEDLIGHTSQUAREINV)*IntegratedPot[1];
ray[0].DRshift-=2.0*IntegratedPot[2]*SPEEDLIGHTSQUAREINV;

ray[0].DTshift-=2.0*IntegratedPot[1]/(SPEEDLIGHT*Hub);

ray[0].D3shift[0]-=2.0*MaxLineOfSight*IntegratedPot[3]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[0]+=2.0*IntegratedPot[4]*SPEEDLIGHTSQUAREINV;

ray[0].D3shift[1]-=2.0*MaxLineOfSight*IntegratedPot[5]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[1]+=2.0*IntegratedPot[6]*SPEEDLIGHTSQUAREINV;

ray[0].D3shift[2]-=2.0*MaxLineOfSight*IntegratedPot[7]*SPEEDLIGHTSQUAREINV;
ray[0].D3shift[2]+=2.0*IntegratedPot[8]*SPEEDLIGHTSQUAREINV;

//0:Total Mag, 1 Converegence, 2 Doppler vel ,3 Total vel with obs, 4=2 convergance for , 5 doppler with velocity shift
ray[0].Magni[0] += -4.*IntegratedPot[0]/(SPEEDLIGHT*SPEEDLIGHT*MaxLineOfSight);
ray[0].Magni[0] += -2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(-ray[0].pot*SPEEDLIGHTSQUAREINV+ray[0].vr/SPEEDLIGHT-2.*IntegratedPot[1]*SPEEDLIGHTSQUAREINV);
ray[0].Magni[0] += 2.*(IntegratedPot[9]-IntegratedPot[10]/MaxLineOfSight)*SPEEDLIGHTSQUAREINV;
ray[0].Magni[1] = 2.*(IntegratedPot[9]-IntegratedPot[10]/MaxLineOfSight)*SPEEDLIGHTSQUAREINV;
ray[0].Magni[2] = -2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(ray[0].vr/SPEEDLIGHT);
ray[0].Magni[3] = ray[0].Magni[0] - 2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(run.Pot_Obs*SPEEDLIGHTSQUAREINV- vr_obs/SPEEDLIGHT) - 2.*vr_obs/SPEEDLIGHT;
ray[0].Magni[4] = ray[0].Magni[2];//Convergance doesn't depend on observer
ray[0].Magni[5] = -2.*(1.-SPEEDLIGHT/(Hub*MaxLineOfSight))*(ray[0].vr-vr_obs)/SPEEDLIGHT;
}
