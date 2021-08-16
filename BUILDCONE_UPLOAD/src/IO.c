#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./vardef.h"

double calc_comoving_distance(double red){
  int i=0;
  double z,result=0.0,dz=red/5000,Om = omega_m,OL = omega_L;
  result=pow(Om*pow(red+1,3)+(red+1)*(red+1)*(1-Om-OL)+OL,-0.5)/2.0+0.5;
  for (z=dz;z<red-dz*0.1;z+=dz,i++) result+=pow(Om*pow(z+1,3)+(z+1)*(z+1)*(1-Om-OL)+OL,-0.5);
  return (c0/hubble_constant)*result*dz; // Mpc/h
  }
float *Read_density_float(char *path,int *Ngrid,float *box_length,float *average_den){
  FILE *save_file;
  int Ngrid3d_lc,num_neg = 0,num_nan = 0;
  float box_lc,dum,av_denc,grid_size;
  float *dens;
  if (!(save_file = fopen (path,"rb"))) {printf("Can't find file,%s\n",path );exit(0);}
  fread(&Ngrid3d_lc,4,1,save_file);
  fread(&box_lc,4,1,save_file);
  fread(&av_denc,4,1,save_file);
  dens = (float *) malloc( Ngrid3d_lc * sizeof(float));
  for (int idx_i = 0; idx_i < Ngrid3d_lc; idx_i++)
  {fread(&dum,4,1,save_file);
    if(isnan(dum)) num_nan++;
    else if (dum<0){num_neg++;}
    dens[idx_i]=  dum;
  }
  fclose(save_file);
  fflush(stdout);
  grid_size = box_lc/ cbrt(Ngrid3d_lc);
  if ((num_nan > 0)||(num_neg > 0)) printf("Warning, there are %d -ve values and %d values,\nI will still run it but you know...\n be carefull bro\n",num_neg,num_nan);
  printf("Density array Read from \t %s\n with grid %d, box size %f and average denisty %f \n",path,(int) cbrt(Ngrid3d_lc), box_lc,av_denc );
  *Ngrid = (int) cbrt(Ngrid3d_lc); *box_length = box_lc;*average_den = av_denc;
  return dens;
}
double *Read_density(char *path,int *Ngrid,float *box_length,float *average_den){
      FILE *save_file;
      int Ngrid3d_lc,num_neg = 0,num_nan = 0;
      float box_lc,dum,av_denc,grid_size;
      double *dens;
      printf("Reading Density ..... \n");fflush(stdout);
      if (!(save_file = fopen (path,"rb"))) {printf("Can't find file,%s\n",path );exit(0);}
      fread(&Ngrid3d_lc,4,1,save_file);
      fread(&box_lc,4,1,save_file);
      fread(&av_denc,4,1,save_file);
      dens = (double *) malloc( Ngrid3d_lc * sizeof(double));
      for (int idx_i = 0; idx_i < Ngrid3d_lc; idx_i++)
        {fread(&dum,4,1,save_file);
          if(isnan(dum)) num_nan++;
          else if (dum<0){num_neg++;}
          dens[idx_i]= (double) dum;
        }
      fclose(save_file);
      fflush(stdout);
      grid_size = box_lc/ cbrt(Ngrid3d_lc);
      if ((num_nan > 0)||(num_neg > 0)) printf("Warning, there are %d -ve values and %d values,\nI will still run it but you know...\n be carefull bro\n",num_neg,num_nan);
      printf("...Done.\nDensity array Read from \t %s\n with grid %d, box size %f and average denisty %f \n",path,(int) cbrt(Ngrid3d_lc), box_lc,av_denc );
      fflush(stdout);
      *Ngrid = (int) cbrt(Ngrid3d_lc); *box_length = box_lc;*average_den = av_denc;
      return dens;
    }

struct io_header_1 load_header_gdt(char *fname,int singlefile){
      /*Reads header from gadget file with path fname"*/
        struct io_header_1 head;
        FILE *fd;
        char blkname[4],buf[200];
        int dummy,SnapFormat=1;
        fflush(stdout);
        if(singlefile == 1)sprintf(buf, "%s", fname);
        else sprintf(buf, "%s.%d", fname, 0);

        if(!(fd = fopen(buf, "r"))){
          printf("can't open header file %s\n", fname);
          exit(0);
        }
        printf("Reading header from `%s' ...", fname);
        fflush(stdout);
        fread(&dummy, sizeof(dummy), 1, fd);
        if (dummy==8) SnapFormat=2; else if (dummy==256) SnapFormat=1; else {printf("Unknown Snapshot format!\nFirst check sum=%i\n",dummy); fclose(fd); exit(0);}
        if (SnapFormat==2){
          fread(blkname,sizeof(char),4,fd);
          fread(&dummy, sizeof(dummy), 1, fd);
          fread(&dummy, sizeof(dummy), 1, fd);
          fread(&dummy, sizeof(dummy), 1, fd);
        }
        fread(&head, sizeof(head), 1, fd);
        fread(&dummy, sizeof(dummy), 1, fd);
        fclose(fd);
        printf(".. done.\n");fflush(stdout);
        return head;
      }
struct survey_functions *Form_surv_struct(char *fname,char delim[]){
          int ncols,nrows,idx=0;
          double chi_max;
          FILE *survey_func_file;
          struct survey_functions *srv_funcs;
          char buf[0x1000];
          Text_file_specs(fname,delim,&ncols,&nrows);
          if(!(srv_funcs = malloc(((nrows+10) * sizeof(struct survey_functions))))) {fprintf(stderr, "failed to allocate memory.\n"); exit(0);}
          if (!(survey_func_file = fopen (fname,"r"))) printf("Can not open survey file %s\n",fname );
          surv_number_pts = nrows;
          while (fgets(buf, sizeof(buf), survey_func_file) != NULL) {
            if ((buf[0] == '#')||(buf[0] == '%')) continue;
            sscanf(buf,"%lf\t%lf\t%lf\t%lf\t%lf",&srv_funcs[idx].redshift,&srv_funcs[idx].b_lin,&srv_funcs[idx].n_z,&srv_funcs[idx].q_mag,&srv_funcs[idx].e_evl);
            idx++;
            }
            fclose(survey_func_file);
            for (int idx_column = 1; idx_column < nrows; idx_column++){
              srv_funcs[idx_column].chi =  calc_comoving_distance(srv_funcs[idx_column].redshift);
              if (srv_funcs[idx_column].b_lin - 1 < 0) chi_neg = srv_funcs[idx_column].chi;//if you're bias goes under zero at some chi the densities are set to zero is that range
              }
            printf("Reading Survey function from %s, the functions are defined up to chi_final %f and redshift  %f\n",fname,srv_funcs[idx-1].chi,srv_funcs[idx-1].redshift);
            fflush(stdout);
            srv_funcs[0].redshift   = 0;
            srv_funcs[0].n_z   = srv_funcs[1].n_z;
            srv_funcs[0].q_mag   = srv_funcs[1].q_mag;
            srv_funcs[0].e_evl   = srv_funcs[1].e_evl;
            srv_funcs[0].b_lin   = srv_funcs[1].b_lin;
            srv_funcs[0].chi   = 0.0;
            if (USESPLINE) Initiate_spline(srv_funcs,nrows);////GSL_SECTIONS
            return srv_funcs;
      }
struct particle_data *Read_Coordinates(char *fname,int n_gal){
      /*Reads co-ordiantes in liger binary files*/
      struct particle_data *par_data;
      FILE * File_r;
      char path[200];
      int nu_nan=0;
      float dist;
      printf("Reading particle coordinates in %s ......",fname);
      fflush(stdout);
      sprintf(path,"%s",fname);
      if(!(par_data = malloc(n_gal * sizeof(struct particle_data)))) printf("Failed to allocate memory\n");
      fflush(stdout);
      if (!(File_r = fopen (path,"r"))) {printf("Can't open file %s to read co-ordinates.\n" ,path);exit(0);}
      fseek(File_r, 8, SEEK_SET);
      for (int idx_read = 0; idx_read < n_gal; idx_read++) {
        fread(&par_data[idx_read].pos[0],4,3,File_r);
        dist = sqrt(par_data[idx_read].pos[0]*par_data[idx_read].pos[0]+par_data[idx_read].pos[1]*par_data[idx_read].pos[1]+par_data[idx_read].pos[2]*par_data[idx_read].pos[2]);
        if(!isfinite(dist)||(dist>24000)){
          nu_nan++;
          par_data[idx_read].pos[0]=0;
          par_data[idx_read].pos[1]=0;
          par_data[idx_read].pos[2]=0;
        }
      }
      printf("...Done \n  Particles have %d NAN values out of %d values\n", nu_nan,n_gal);
      fclose(File_r);
      fflush(stdout);
      return par_data;
    }
struct magni *Read_magnification(char *path, int n_gal,int r_max){
      FILE * File_mag;
      char magni_path[200];
      int ngal_c,n_neg=0,n_1plus=0;
      float rmax_c;
      struct magni *mag;
      sprintf(magni_path,"%s",path);
      if (!(File_mag = fopen (magni_path,"r"))) {printf("Can't open Mag file%s\n",magni_path );exit(0);}
      printf("Reading Magnification ...\n");
      fflush(stdout);
      fread(&rmax_c,4,1,File_mag);
      fread(&ngal_c,4,1,File_mag);
      if ((rmax_c==r_max)&&(ngal_c==n_gal)) printf("Mag Files Have  Correct Header\n");else printf("Mag FIle header doesn't match to coordinate file headers\n");
      if(!(mag = malloc(n_gal * sizeof(struct magni)))) printf("Failed to allocate memory\n");
      fflush(stdout);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) fread(&mag[idx_r].mag_tot,4,1,File_mag);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) fread(&mag[idx_r].mag_conv,4,1,File_mag);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) fread(&mag[idx_r].mag_doppler,4,1,File_mag);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) fread(&mag[idx_r].mag_tot_obs,4,1,File_mag);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) fread(&mag[idx_r].mag_doppler_obs,4,1,File_mag);
      printf("...Done\t\n");
      fflush(stdout);
      fclose(File_mag);
      fflush(stdout);
      return mag;
    }
struct Time_data *Read_dlna(char *path, int n_gal,int r_max){
      FILE * File_TIME;
      char time_path[200];
      int ngal_c;
      float rmax_c;
      struct Time_data *time_data,*buffer;
      sprintf(time_path,"%s",path);
      if (!(File_TIME = fopen (time_path,"r"))) {printf("Can't open Time file%s\n",time_path );exit(0);}
      printf("Reading Time_data ...\n");
      fflush(stdout);
      fread(&rmax_c,4,1,File_TIME);
      fread(&ngal_c,4,1,File_TIME);
      if ((rmax_c==r_max)&&(ngal_c==n_gal)) printf("Time_data Files Have  Correct Header\n");else printf("Time Data FIle header doesn't match to coordinate file headers\n");
      if(!(time_data = malloc(n_gal * sizeof(struct Time_data)))) printf("Failed to allocate memory\n");
      if(!(buffer = malloc(2 * sizeof(struct Time_data)))) printf("Failed to allocate memory\n");//We don't read the time distance because we don't need it
      fflush(stdout);
      for (int idx_r = 0; idx_r < n_gal; idx_r++) {
        fread(&buffer[0].dlna[0],4,4,File_TIME);
        fread(&time_data[idx_r].dlna[0],4,4,File_TIME);
      }
      int rand_int;
      printf("DLNA CHECK\n");
      for (int idx_r = 0; idx_r < 10; idx_r++) {
        rand_int = rand() % n_gal;
        for (int idx_print = 0; idx_print < 4; idx_print++)   printf("%e\t",time_data[rand_int].dlna[idx_print]);
        printf(" %f\n",time_data[rand_int].dlna[1]/time_data[rand_int].dlna[3]);
      }
      printf("...Done\t\n");
      fflush(stdout);
      fclose(File_TIME);
      fflush(stdout);
      return time_data;
    }
struct dens_com *Read_nbar_text(char *save_path,int *array_size){
      struct dens_com *n_bar;
      FILE *f = fopen(save_path, "r");if (f == NULL) {printf("Error opening file!\n");exit(0);} else printf("File opened for reading density\n" );
      int nu_of_pts;
      fscanf(f, "%d", &nu_of_pts);
      if (!(n_bar =  malloc( nu_of_pts * sizeof(struct dens_com)) )) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
      // fprintf(f, "%f\n", grid_size);
      for (int idx_w = 0; idx_w < nu_of_pts; idx_w++) {fscanf(f, "%le\t%le", &n_bar[idx_w].chi,&n_bar[idx_w].n_bar);
      // if (depugging) printf("%le\t%le\n", n_bar[idx_w].chi,n_bar[idx_w].n_bar);
        }

      fclose(f);
      fflush(stdout);
      *array_size = nu_of_pts;
      printf("n_bar Text array read at \t %s\n",save_path );fflush(stdout);
      return n_bar;
    }

void File_Specs(char *particle_path,int *n_gal,float *r_max){
  /*Reads header in liger binary files
  suffix:set extension
  particle_path: directory of binary file
  n_gal and r_max containers for the header
  */
  printf("Reading LIGER file header....");
  fflush(stdout);
  int n_gal_f_r,file_size_nu;
  long file_size_bytes;
  float r_max_f_r;
  FILE * pFile;
  char path[200];
  fflush(stdout);
  sprintf(path,"%s",particle_path);
  fflush(stdout);
  //Reading the full_coordinates file
  if (!(pFile = fopen (path,"r"))) {printf("Can't read specs %s\n",path );exit(0);}
  fseek(pFile, 0L, SEEK_END);
  file_size_bytes = ftell(pFile);
  file_size_nu = (file_size_bytes)/4 - 2;
  rewind(pFile);
  fread(&r_max_f_r,4,1,pFile);
  fread(&n_gal_f_r,4,1,pFile);
  if (n_gal_f_r == file_size_nu/3) {
    printf("...Done\n");
    printf ("LIGTCONE Size: %f and NUMBER of particles: %d \n",r_max_f_r,n_gal_f_r);
    fflush(stdout);
  } else {
    printf("File size doesn't match number of particles in headwer,Corrupt file,Number = %d\t Number in files %d \t\n",n_gal_f_r,file_size_nu/3);
  }
  fclose(pFile);
  //variables for output
  *n_gal=n_gal_f_r;
  *r_max=r_max_f_r;
}
void Text_file_specs(char *fname,char delim[],int *ncol,int *nrow){
  int nrow_dum = 0,ncol_dum = 0;
  FILE *fp;
  char c;
  char *result = NULL;
  if (!(fp = fopen (fname,"r"))) printf("Can not open text file %s\n",fname );
  for (c = getc(fp); c != EOF; c = getc(fp)){

   if ((c ==delim)&&(nrow_dum==0)) ncol_dum++;
   if (c =='\n') nrow_dum++;
    }
  // printf("File %s has %d rows and %d colums  \n",fname,nrow_dum,ncol_dum );
  *ncol= ncol_dum+1;
  *nrow= nrow_dum;
}
void Initiate_spline(struct survey_functions *surv,int nrows){
  int n_points=0,idx_spline_st=0;
  printf("Initiating Cubic Splines with %d elements ...\t",nrows);
  if (!(surv_splines.CHI_interp=(double*) malloc((nrows-idx_spline_st)*sizeof(double)))){	printf("\nAllocation of memory failed.\n");exit(0);}
  if (!(surv_splines.Func_CHI=(double*) malloc((nrows-idx_spline_st)*sizeof(double)))){	printf("\nAllocation of memory failed.\n");exit(0);}

  surv_splines.acc_red = gsl_interp_accel_alloc();
  surv_splines.acc_b_lin = gsl_interp_accel_alloc();
  surv_splines.acc_n_z = gsl_interp_accel_alloc();
  surv_splines.acc_q_mag = gsl_interp_accel_alloc();
  surv_splines.acc_e_bias = gsl_interp_accel_alloc();
  //allocate spine
  surv_splines.spline_red = gsl_spline_alloc (gsl_interp_cspline, nrows-idx_spline_st);
  surv_splines.spline_b_lin = gsl_spline_alloc (gsl_interp_cspline, nrows-idx_spline_st);
  surv_splines.spline_n_z = gsl_spline_alloc (gsl_interp_cspline, nrows-idx_spline_st);
  surv_splines.spline_q_mag = gsl_spline_alloc (gsl_interp_cspline, nrows-idx_spline_st);
  surv_splines.spline_e_bias = gsl_spline_alloc (gsl_interp_cspline, nrows-idx_spline_st);
  //redshift
  {  for (int  idx_bar = 0; idx_bar < nrows-idx_spline_st; idx_bar++) {
      surv_splines.CHI_interp[idx_bar] = surv[idx_bar+idx_spline_st].chi;
      surv_splines.Func_CHI[idx_bar] = surv[idx_bar+idx_spline_st].redshift;
      if ((surv_splines.CHI_interp[idx_bar]==0.0)&&(idx_bar>2))surv_splines.CHI_interp[idx_bar] = 2*surv_splines.CHI_interp[idx_bar-1] - surv_splines.CHI_interp[idx_bar-2];
    }
    gsl_spline_init (surv_splines.spline_red, surv_splines.CHI_interp, surv_splines.Func_CHI, nrows-idx_spline_st);
  }
  {  for (int  idx_bar = 0; idx_bar < nrows-idx_spline_st; idx_bar++) {
      surv_splines.CHI_interp[idx_bar] = surv[idx_bar+idx_spline_st].chi;
      surv_splines.Func_CHI[idx_bar] = surv[idx_bar+idx_spline_st].b_lin;
      if ((surv_splines.CHI_interp[idx_bar]==0.0)&&(idx_bar>2))surv_splines.CHI_interp[idx_bar] = 2*surv_splines.CHI_interp[idx_bar-1] - surv_splines.CHI_interp[idx_bar-2];
    }
    gsl_spline_init (surv_splines.spline_b_lin, surv_splines.CHI_interp, surv_splines.Func_CHI, nrows-idx_spline_st);
  }
  {  for (int  idx_bar = 0; idx_bar < nrows-idx_spline_st; idx_bar++) {
      surv_splines.CHI_interp[idx_bar] = surv[idx_bar+idx_spline_st].chi;
      surv_splines.Func_CHI[idx_bar] = surv[idx_bar+idx_spline_st].n_z;
      if ((surv_splines.CHI_interp[idx_bar]==0.0)&&(idx_bar>2))surv_splines.CHI_interp[idx_bar] = 2*surv_splines.CHI_interp[idx_bar-1] - surv_splines.CHI_interp[idx_bar-2];
    }
    gsl_spline_init (surv_splines.spline_n_z, surv_splines.CHI_interp, surv_splines.Func_CHI, nrows-idx_spline_st);
  }
  {  for (int  idx_bar = 0; idx_bar < nrows-idx_spline_st; idx_bar++) {
      surv_splines.CHI_interp[idx_bar] = surv[idx_bar+idx_spline_st].chi;
      surv_splines.Func_CHI[idx_bar] = surv[idx_bar+idx_spline_st].q_mag;
      if ((surv_splines.CHI_interp[idx_bar]==0.0)&&(idx_bar>2))surv_splines.CHI_interp[idx_bar] = 2*surv_splines.CHI_interp[idx_bar-1] - surv_splines.CHI_interp[idx_bar-2];
    }
    gsl_spline_init (surv_splines.spline_q_mag, surv_splines.CHI_interp, surv_splines.Func_CHI, nrows-idx_spline_st);
  }
  //Evloitionary bias
  {  for (int  idx_bar = 0; idx_bar < nrows-idx_spline_st; idx_bar++) {
      surv_splines.CHI_interp[idx_bar] = surv[idx_bar+idx_spline_st].chi;
      surv_splines.Func_CHI[idx_bar] = surv[idx_bar+idx_spline_st].e_evl;
      if ((surv_splines.CHI_interp[idx_bar]==0.0)&&(idx_bar>2))surv_splines.CHI_interp[idx_bar] = 2*surv_splines.CHI_interp[idx_bar-1] - surv_splines.CHI_interp[idx_bar-2];
    }
    gsl_spline_init (surv_splines.spline_e_bias, surv_splines.CHI_interp, surv_splines.Func_CHI, nrows-idx_spline_st);
  }
  free(surv_splines.CHI_interp);
  free(surv_splines.Func_CHI);
  printf("..Done\n");
}
void Free_spline(){
  gsl_interp_accel_free(surv_splines.acc_red);
  gsl_interp_accel_free(surv_splines.acc_b_lin);
  gsl_interp_accel_free(surv_splines.acc_n_z);
  gsl_interp_accel_free(surv_splines.acc_q_mag);
  gsl_interp_accel_free(surv_splines.acc_e_bias);
  //allocate spine
  gsl_spline_free(surv_splines.spline_red );
  gsl_spline_free(surv_splines.spline_b_lin );
  gsl_spline_free(surv_splines.spline_n_z );
  gsl_spline_free(surv_splines.spline_q_mag );
  gsl_spline_free(surv_splines.spline_e_bias);

}
void Save_array(char *save_path,float *array, int array_size,float box_len,float av_den){
      FILE *save_file;
      printf("Saving Binary Array.....");
      fflush(stdout);
      if(!(save_file = fopen (save_path,"wb"))) {printf("Can't save file%s\n",save_path);exit(0);}
      fwrite(&array_size,4,1,save_file);
      fwrite(&box_len,4,1,save_file);
      fwrite(&av_den,4,1,save_file);
      for (int idx_i = 0; idx_i < array_size; idx_i++) fwrite(&array[idx_i],4,1,save_file);
      fclose(save_file);
      printf("...Done\n Saved at \t %s.\n",save_path );
      fflush(stdout);
    }
void Save_array_text(char *save_path,double *array, int array_size,float grid_size,double *k_bin){
      FILE *f = fopen(save_path, "w");if (f == NULL) {printf("Error opening file!\n");exit(0);} else printf("File opened for saving text array\n" );
      fprintf(f, "%d\t%f\n", array_size,grid_size);
      // fprintf(f, "%f\n", grid_size);
      for (int idx_w = 0; idx_w < array_size; idx_w++) fprintf(f, "%e\t%e\n", array[idx_w],k_bin[idx_w]);
      fclose(f);
      fflush(stdout);
      printf("Text array saved at \t %s\n",save_path );fflush(stdout);
    }
void Save_nbar_text(char *save_path,struct dens_com *n_bar, int array_size){
      FILE *f = fopen(save_path, "w");if (f == NULL) {printf("Error opening file!\n");exit(0);} else printf("File opened for writing density\n" );
      fprintf(f, "%d\t\n", array_size);
      // fprintf(f, "%f\n", grid_size);
      for (int idx_w = 0; idx_w < array_size; idx_w++) fprintf(f, "%le\t%le\n", n_bar[idx_w].chi,n_bar[idx_w].n_bar);
      fclose(f);
      fflush(stdout);
      printf("Text array saved at \t %s\n",save_path );fflush(stdout);
    }
void Read_simulation_parameters(char *music_base){
  /*Wraper around header reader to get the density parameters*/
  char *header_path[200];
  struct io_header_1 main_head;
  main_head = load_header_gdt(music_base,1);
  omega_L = main_head.OmegaLambda;
  omega_m = main_head.Omega0;
  n_sim_bar = main_head.npartTotal[1]/(main_head.BoxSize*main_head.BoxSize*main_head.BoxSize);
  printf("Simulation Parameters (Omega_m , Omega_L, n_sim_bar)\n\t\t(%f, %f, %f)\n",omega_m, omega_L, n_sim_bar);
}

// struct io_header_1 header_spare(int ngal,)
