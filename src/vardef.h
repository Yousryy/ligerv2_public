/*  This file is part of Liger, see main.c and the User guide.
    Copyright (C) 2017  Mikolaj Borzyszkowski, Argelander Institut für Astronomie, University of Bonn
    Copyright (C) 2021  Mohamed Elkhashab, Argelander Institut für Astronomie, University of Bonn
*/

#ifndef _VARDEF_H
#define _VARDEF_H

//This struct holds the main parameters of the simulation.
struct header_snapshot
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

//This struct holds the parameters form the Parameterfile.
struct header_runpara
{
  char outname[500];
  char outpath[500];
  char gravpath[500];
  char simpath[500];
  char simpath_gal[500];
  int NSnap;
  int ifOutputPotential;
  int Nthreads;
  int Ngrid;
  int NPartInMem;
  double maxdist;
  double mindist;
  double ObsPos[3];
  double OutBuf;
  //Internally found variables
  double ObsVel[3];
  double Box_size;
  double Pot_Obs;
} run;

//This struct is used for the distribution of matter particles and galaxies.
struct particle_data
{
  float Pos[3];
  float Vel[3];
};

//This struct summarizes the final result of Liger.
struct light_cone_object
{
  float SimPos[3];
  float VelPos[3];
  float VelPos_obs[3];
  float RelPos[3];
  float RelPos_obs[3];
  float Magni[6];
  float TimePosition[8];// First four are time positions in distance unit , second 4 are logarthmic lna

//  int id;
};


//This struct is used to define the line-of-sight which should be integrated and holds the final corrdiante transfomation.
struct light_rays_path {
  double start[3];
  double end[3];
  float normvec[3];
  int CellsAlongPath;
  int *PathCell;
  float *PathIntersec;
  float PathIntersec_max_line;
  float dist;
  float pot;
  float TimeDistance;
  float dlna;
  float vr;
  float vrHub;
  float D3shift[3];
  float DRshift;
  float DTshift;
  float Magni[6];
};

//Here we store the potantial together with its derivatives.
struct potential_container
{
  float val;
  float timeD;
  float firD[3];
  float secD[3];
  float croD[3];
};
//Global Variables
double *rangemin,*rangemax,*Hubarray;
int SNAP_ORDER;
#ifdef MEASURETIME
  time_t TotalTime;
  time_t PartTime[2];
  time_t IOTime[2];
#endif
//long Ended_before_crossing_light_cone=0;
#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })
#define MAXPATHLENGTH 200 //Maximum total character lenght of the path to any file.
#define CONSPI (3.141592654)
#define SPEEDLIGHT (2.99724e5)
#define SPEEDLIGHTSQUAREINV (1.11265e-11)

#endif
