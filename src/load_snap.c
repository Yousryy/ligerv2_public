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

//  This file contains the routines which read the simulation data and the list of galaxies.
//  These routines only redirect their input, this is doe to facilitate the use of custom routines to handle differet format of input data.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "./vardef.h"
#include "./load_snap_gdt.h"

struct header_snapshot load_snapshot_header(int FileNr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);	fflush(stdout);
return load_header_gdt(fname);
}

struct header_snapshot load_galaxies_header(int FileNr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);	fflush(stdout);
return load_header_gdt(fname);
}

int load_snapshot(int FileNr, int files,struct particle_data **Pptr,int read_velocity){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath,FileNr);	fflush(stdout);
#ifdef MEASURETIME
 IOTime[0]-=time(NULL);
 int result=load_snapshot_gdt(fname,files,Pptr,read_velocity);
 IOTime[0]+=time(NULL);
 return result;
#else
 return load_snapshot_gdt(fname,files,Pptr,read_velocity);
#endif
}

//Currently for Gadget only files=1 is possible.
int load_galaxies_partial(int FileNr,int files,int NStart,int NSelect,struct particle_data **Pptr){
char fname[MAXPATHLENGTH];
sprintf(fname,"%s%03d",run.simpath_gal,FileNr);	fflush(stdout);
#ifdef MEASURETIME
 #pragma omp atomic
  IOTime[1]-=time(NULL);
 int result=load_snapshot_gdt_partial(fname,files,NStart,NSelect,Pptr);
 #pragma omp atomic
  IOTime[1]+=time(NULL);
 return result;
#else
 return load_snapshot_gdt_partial(fname,files,NStart,NSelect,Pptr);
#endif
}
int find_snapshot_ordering(){
  //Check how the snapshots are ordered in redshift to interpolate between snapshots in the correct "Direction"
  struct header_snapshot Header_i,Header_f;
  float z_i,z_f;
  Header_i=load_snapshot_header(0);
  Header_f=load_snapshot_header(run.NSnap-1);
  z_i = Header_i.redshift;
  z_f = Header_f.redshift;
  if (z_f>z_i){
    printf("Input is ordered in acsending order in z  %03d to %03d = %.1f to %.1f \n",0,run.NSnap-1,z_i,z_f);
    fflush(stdout);
    return 0;
  }else if (z_i>z_f){
    printf("Input is ordered in decsending order in z  %03d to %03d = %.1f to %.1f \n",0,run.NSnap-1,z_i,z_f);
    fflush(stdout);
    return 1;
  }else{
    printf("First and last file have same redshift No evloution  %03d to %03d = %.1f to %.1f \n",0,run.NSnap-1,z_i,z_f);
    exit(0);
  }
}
