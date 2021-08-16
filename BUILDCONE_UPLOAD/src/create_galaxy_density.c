// Add option to save individual things and plot

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "./vardef.h"
#include "./IO.h"
#include "./create_galaxy_density.h"





void Build_Cone(char *base_liger,char *base_music,char *surv_path,char *save_path,int Ngrid){
  /*Builds light cone desnity according to liger paper
  mode specifies kind of effect:
  0:Realspace
  1:GRRSD
  2:GRRSD Observer
  3:GRRSD With Constant Q in vardef
  4:GRRSD With Constant Q = 0
  5:vRSD
  6:vRSD Observor
  7:GRRSD Observer with Constant Q in vardef
  8:GRRSD with Constant Q =1
  9:GRRSD Observer with Constant Q =1
  10:Produces all densities in list defined in vardef
  Input are directories passed by main given in command line
  */
  float Q_const = run.Q_const;
  char redshift_path_suffix[10][20] = {"_realspace.dat", "_GRRSD_obs.dat","_GRRSD_obs.dat"  ,"_GRRSD_obs.dat" ,"_GRRSD.dat", "_GRRSD.dat"    ,  "_GRRSD.dat" , "_GRRSD.dat" ,  "_vRSD.dat" , "_vRSD_obs.dat"};
  char save_path_suffix[10][20]     = {"_realspace.dat", "_GRRSD_obs.dat","_GRobsQconst.dat","_GRobsQ1.dat"   ,"_GRRSD.dat", "_GRQconst.dat" ,  "_GRQ1.dat"  , "_GRQ0.dat " ,  "_vRSD.dat" , "_vRSD_obs.dat"};
  sprintf(save_path_suffix[2],"%s%d%s","_GRobsQ",(int) Q_const,".dat");//Change _GRobsQconst with const value
  sprintf(save_path_suffix[5],"%s%d%s","_GRQ",(int) Q_const,".dat");//Change GRQ with const value
  int mode_mag_arr[10]                  = {0  , 3,  3,  3,  2,  2, 2,  0, 0, 0};
  int mode_evol_arr[10]                 = {0  , 4,  4,  4,  5,  5, 5 , 5, 7, 6};
  float q_values_arr[10]                = {0. ,-1., 5., 1.,-1., 5.,1., 0, 0, 0};
  q_values_arr[2] = Q_const;
  q_values_arr[5] = Q_const;
  char surv_func_path[200],real_path[200],magni_path[200],time_path[200],redshift_path[200],test_path[200];
  int n_gal_r,n_gal_s,Ngrid3d,mode_mag,mode_evol;
  float r_max_lgr,l_box, gridsize, chi_max_par,av_den,q_value;
  int rand_int;//conssitency
  float r_zbar, r_z;//consistency
  int wanted_dens[run.nu_of_cones];
  for (int idx_d = 0; idx_d < run.nu_of_cones; idx_d++) wanted_dens[idx_d] = run.dens_to_compute[idx_d];
  struct survey_functions *survey;
  struct particle_data *real_par,*shift_par;
  struct magni *mag;
  struct Time_data *T_dlna;
  struct io_header_1 main_head;
  Ngrid3d = Ngrid*Ngrid*Ngrid;
  float *gal_density;
  sprintf(magni_path,"%s%s", base_liger,"_magni.dat");
  sprintf(time_path,"%s%s", base_liger,"_timedata.dat");
  sprintf(real_path,"%s%s", base_liger,"_realspace.dat");
  sprintf(surv_func_path, "%s",surv_path);
  float factor_box_size = run.factor_box_size;
  //Reading Cosmo parameters
  // main_head = load_header_gdt(base_music,1);
  // omega_L = main_head.OmegaLambda;
  // omega_m = main_head.Omega0;
  // n_sim_bar = main_head.npartTotal[1]/(main_head.BoxSize*main_head.BoxSize*main_head.BoxSize);
  //Checking that no files are overwritten
  int exist_test = 1,idx_check;
  printf("Checking that no files will be overwritten.\n");
  {
  for (int counter = 0;counter<sizeof(wanted_dens)/sizeof(wanted_dens[0]) ; counter++){
    idx_check  = wanted_dens[counter];
    sprintf(test_path,"%s%s_%d%s", save_path,"dens", Ngrid,save_path_suffix[idx_check]);
    if( access( test_path, F_OK ) != -1 ){
    printf("%s Already built, will skip this cone \n",test_path);fflush(stdout);
    }else{
      printf("%s Will be built , \n",test_path);fflush(stdout);
      exist_test = 0;
    }
  }
        if (exist_test == 1) {
      printf("ALL cones requested already build, please delete them to write\n");fflush(stdout);
      return;
    }
  }
  printf("Simulation Parameters used for distance computations(Omega_m , Omega_L, n_sim_bar)\n\t\t(%f, %f, %f)\n",omega_m, omega_L, n_sim_bar);
  survey = Form_surv_struct(surv_func_path,' ');
  int suf_counter = 0 , mode_counter = 0 , idx_steps=0 , dens_build_nu=0 , idx_next , idx_bef;
  float  *real_dens,*redshift_dens,*E_bias_dens,*Q_bias_dens;
  File_Specs(real_path,&n_gal_r,&r_max_lgr);
  mag = Read_magnification(magni_path,n_gal_r,r_max_lgr);
  real_par = Read_Coordinates(real_path,n_gal_r);
  T_dlna = Read_dlna(time_path,n_gal_r,r_max_lgr);
  l_box = factor_box_size*r_max_lgr;
  gridsize = (l_box / (float) Ngrid);
  l_box += gridsize;
  gridsize = (l_box / (float) Ngrid);
  real_dens = CIC(real_par,mag,T_dlna,survey,gridsize,Ngrid,n_gal_r,0);

  dens_build_nu = sizeof(wanted_dens)/sizeof(wanted_dens[0]);
  for (int counter = 0;counter< dens_build_nu ; counter++){
    suf_counter  = wanted_dens[counter];
    if (suf_counter == 0) continue;
    mode_mag = mode_mag_arr[suf_counter];
    mode_evol = mode_evol_arr[suf_counter];
    q_value = q_values_arr[suf_counter];
    if (counter+1<dens_build_nu) idx_next = counter+1;
    else idx_next = 0;
    if (counter-1 > 0) idx_bef = counter-1;
    else idx_bef = 0;
    sprintf(save_path,"%s%s_%d%s", save_path,"dens", Ngrid,save_path_suffix[suf_counter]);
    sprintf(redshift_path,"%s%s", base_liger,redshift_path_suffix[suf_counter]);
    if( access( save_path, F_OK ) != -1 ){
      printf("%s Already built, will skip this cone \n",save_path);fflush(stdout);
      continue;
    }

    // if ((idx_steps>0)&&(strcmp(redshift_path_suffix[suf_counter],redshift_path_suffix[wanted_dens[idx_bef]])))
    if ((idx_steps==0)||(mode_evol != mode_evol_arr[wanted_dens[idx_bef]]))
    { printf("Recreating Densities %i %i %i",idx_bef , mode_evol,mode_evol_arr[wanted_dens[idx_bef]]);
    File_Specs(redshift_path,&n_gal_s,&r_max_lgr);
    if (n_gal_s != n_gal_r){printf("Warning:Incompatable number of particles between files, Please make sure all files were produced in same run\n");}
    shift_par = Read_Coordinates(redshift_path,n_gal_s);
    chi_max_par = Max_av_shift(n_gal_r,real_par,shift_par);
    redshift_dens = CIC(shift_par,mag,T_dlna,survey,gridsize,Ngrid,n_gal_r,1);
    E_bias_dens =  CIC(shift_par,mag,T_dlna,survey,gridsize,Ngrid,n_gal_r,mode_evol);
  }
  av_den = (float)n_gal_r / (4/3*PI*chi_max_par*chi_max_par*chi_max_par);
  if (mode_mag>0 ) {
    if ((idx_steps==0)||(mode_mag != mode_mag_arr[wanted_dens[idx_bef]])) {
      Q_bias_dens = CIC(shift_par,mag,T_dlna,survey,gridsize,Ngrid,n_gal_r,mode_mag);
    }
    gal_density = Create_galaxy_density(real_dens,redshift_dens,E_bias_dens,Q_bias_dens,survey,Ngrid,gridsize,chi_max_par,suf_counter,q_value);
  }else{
    gal_density = Create_galaxy_density(real_dens,redshift_dens,E_bias_dens,E_bias_dens,survey,Ngrid,gridsize,chi_max_par,suf_counter,q_value);
  }
  if (run.SN) ADD_shot(gal_density, Ngrid, l_box);
  Radial_bins(gal_density,survey,0.9, 1.8,gridsize,Ngrid);
  Save_array(save_path,gal_density,Ngrid3d,l_box,av_den);
  free(gal_density);
  if ((idx_next==0) || (mode_evol != mode_evol_arr[wanted_dens[idx_next]])) {
    free(redshift_dens);free(shift_par);free(E_bias_dens);
    n_gal_s = 0;chi_max_par = 0;}
    if ((idx_next==0) || (mode_mag != mode_mag_arr[wanted_dens[idx_next]])) {
      if(mode_mag>0) free(Q_bias_dens);
    }

    idx_steps+=1;
  }
  printf("All Shifted Cones Built... Yep Yep\n");
  printf("Hurrayyyyyy.\n");
  sprintf(save_path,"%s%s_%d%s", save_path,"dens", Ngrid,save_path_suffix[0]);
  fflush(stdout);
  if( access( save_path, F_OK ) != -1 ){
    printf("%s Already built, will skip this cone \n",save_path);fflush(stdout);
  }else{
    if (wanted_dens[0] == 0) {
      mode_mag = mode_mag_arr[0];
      mode_evol = mode_evol_arr[0];
      q_value = q_values_arr[0];
      chi_max_par = Max_av_shift(n_gal_r,real_par,real_par);
      gal_density = Create_galaxy_density(real_dens,real_dens,real_dens,real_dens,survey,Ngrid,gridsize,chi_max_par,0,q_value);
      if (run.SN) ADD_shot(gal_density, Ngrid, l_box);
      Radial_bins(gal_density,survey,0.9, 1.8,gridsize,Ngrid);
      Save_array(save_path,gal_density,Ngrid3d,l_box,av_den);
      free(gal_density);
    }
  }
    free(real_par);
    free(real_dens);
  if (USESPLINE) Free_spline();
  // if (only_dens==1) {printf("Only galaxy densities requested, will exit\n");exit(0);}


}
float Max_av_shift(int n_gal,struct particle_data *real_par,struct particle_data *shift_par ){
      /* Compute max and average difference between two cartiesean arrays and returns maximum value of array*/
      float max=0.0,avg=0.0,chi_particle_r,chi_particle_s,chi_max=0,chi_shift;
      printf("Computing Maximum Radius and Maximum Shift......");
      fflush(stdout);
      for (int idx_particle = 0; idx_particle < n_gal; idx_particle++) {
        chi_particle_r = sqrt(real_par[idx_particle].pos[0]*real_par[idx_particle].pos[0] + real_par[idx_particle].pos[1]*real_par[idx_particle].pos[1] + real_par[idx_particle].pos[2]*real_par[idx_particle].pos[2]);
        chi_particle_s = sqrt(shift_par[idx_particle].pos[0]*shift_par[idx_particle].pos[0] + shift_par[idx_particle].pos[1]*shift_par[idx_particle].pos[1] + shift_par[idx_particle].pos[2]*shift_par[idx_particle].pos[2]);
        chi_shift = fabs(chi_particle_r-chi_particle_s);
        if(chi_particle_r > chi_max) chi_max = chi_particle_r;
        if(chi_particle_s > chi_max) chi_max = chi_particle_s;
        if(isnan(chi_shift)) {printf("problem here%d,%f,%f,%f,%f,%f,%f\n",idx_particle,chi_particle_r,chi_particle_s,chi_shift,shift_par[idx_particle].pos[0], shift_par[idx_particle].pos[1], shift_par[idx_particle].pos[2] );fflush(stdout);}
        avg+=chi_shift;
        if (chi_shift>max) max=chi_shift;
      }
      // printf("Max Shift,AvgShift%f\t%f\t for %d\n",max,avg,n_gal );
      fflush(stdout);
      avg/=(float) n_gal;
      printf("....Done.\nMax Shift,AvgShift with maximum radius \n%f\t%f\t%f for %d\n",max,avg,chi_max,n_gal);
      fflush(stdout);
      return chi_max;
}
float *CIC(struct particle_data *par_struct,struct magni *magn_arr,struct Time_data *T_DLNA,struct survey_functions *srv_funcs,float grid_size,int Ngrid, int n_gal,int mode){
      //r_or_s =0 for weighted real_space density,1 for weighted shifted density ,2 for bias weighted real shifted ,3 for unweighted real density
      /*Cloud in cell algorithim for liger densities and general particle distributions
        Weights are defined in printfs modes 2 and 3 are deprecated*/
      if (mode == 0) printf("Cloud in Cell for unweighted number density for real space .......");
      else if (mode == 1) printf("Cloud in Cell for unweighted number density for redshift space   .......");
      else if (mode == 2) printf("Cloud in Cell with weight = M_GRRSD(z)......." );
      else if (mode == 3) printf("Cloud in Cell with weight = M_GRRSD_OBS(z)......." );
      else if (mode == 4) printf("Cloud in Cell with weight = dlna_GRobs(z)......." );
      else if (mode == 5) printf("Cloud in Cell with weight = dlna_GR(z)......." );
      else if (mode == 6) printf("Cloud in Cell with weight = dlna_vobs(z)......." );
      else if (mode == 7) printf("Cloud in Cell with weight = dlna_v(z)......." );
      else {printf("Unallowed choice 0<=mode<=7, Due to your greed, the code will perish\n");exit(0);}
      fflush(stdout);
      int idx_i,idx_j,idx_k,idx_i_cl,idx_j_cl,idx_k_cl, Ngrid3d = Ngrid*Ngrid*Ngrid;
      double d_x,d_y,d_z,x_c,y_c,z_c,x_p,y_p,z_p,x_p_w,y_p_w,z_p_w,particle_weight=1.0,n_g_part,some_r,chi_par,chi_par_w;
      float Box = (float) Ngrid*grid_size ;
      struct particle_data *main_struct,*weight_struct;
      float *density;
      //testing Variables
      double nu_par=0.0,one_par = 0.0;
      main_struct = par_struct;
      if (!(density = (float *) calloc( Ngrid3d , sizeof(float)) )) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
      for (int idx_f = 0; idx_f < n_gal; idx_f++) {
        x_p = main_struct[idx_f].pos[0];
        y_p = main_struct[idx_f].pos[1];
        z_p = main_struct[idx_f].pos[2];
        if      (mode==2)  particle_weight = magn_arr[idx_f].mag_tot;
        else if (mode==3)  particle_weight = magn_arr[idx_f].mag_tot_obs;
        else if (mode==4)  particle_weight = T_DLNA[idx_f].dlna[3];
        else if (mode==5)  particle_weight = T_DLNA[idx_f].dlna[2];
        else if (mode==6)  particle_weight = T_DLNA[idx_f].dlna[1];
        else if (mode==7)  particle_weight = T_DLNA[idx_f].dlna[0];
        else particle_weight = 1.0;
        chi_par = sqrt(x_p*x_p + y_p*y_p + z_p*z_p);
        x_p += Box/2;
        y_p += Box/2;
        z_p += Box/2;
        idx_i =x_p/grid_size;
        idx_j =y_p/grid_size;
        idx_k =z_p/grid_size;
        x_c =(double) (idx_i + 1 )*grid_size;
        y_c =(double) (idx_j + 1 )*grid_size;
        z_c =(double) (idx_k + 1 )*grid_size;
        d_x = (x_c - x_p)/grid_size;
        d_y = (y_c - y_p)/grid_size;
        d_z = (z_c - z_p)/grid_size;
        if (d_x < 0.5) {
          d_x+=0.5;
          idx_i_cl = idx_i + 1;
          if (idx_i == Ngrid-1) idx_i_cl = 0;
        }else{
          d_x = 1.5 - d_x;
          idx_i_cl = idx_i - 1;
          if (idx_i == 0) idx_i_cl = Ngrid -1;
        }
        if (d_y < 0.5){
          d_y+=0.5;
          idx_j_cl = idx_j + 1;
          if (idx_j == Ngrid-1) idx_j_cl = 0;
        }else{
          d_y = 1.5 - d_y;
          idx_j_cl = idx_j - 1;
          if (idx_j == 0) idx_j_cl = Ngrid -1;
        }
        if (d_z < 0.5){
          d_z+=0.5;
          idx_k_cl = idx_k + 1;
          if (idx_k == Ngrid-1) idx_k_cl = 0;
        } else {
          d_z = 1.5 - d_z;
          idx_k_cl = idx_k - 1;
          if (idx_k == 0) idx_k_cl = Ngrid -1;
        }
        density[Shift_koor(idx_i      ,idx_j      ,idx_k      ,Ngrid,0)] +=  particle_weight*d_x    *d_y    *d_z      ;
        density[Shift_koor(idx_i_cl   ,idx_j      ,idx_k      ,Ngrid,0)] +=  particle_weight*(1-d_x)*d_y    *d_z      ;
        density[Shift_koor(idx_i      ,idx_j_cl   ,idx_k      ,Ngrid,0)] +=  particle_weight*d_x    *(1-d_y)*d_z      ;
        density[Shift_koor(idx_i      ,idx_j      ,idx_k_cl   ,Ngrid,0)] +=  particle_weight*d_x    *d_y    *(1-d_z)  ;
        density[Shift_koor(idx_i_cl   ,idx_j_cl   ,idx_k      ,Ngrid,0)] +=  particle_weight*(1-d_x)*(1-d_y)*d_z      ;
        density[Shift_koor(idx_i_cl   ,idx_j      ,idx_k_cl   ,Ngrid,0)] +=  particle_weight*(1-d_x)*d_y    *(1-d_z)  ;
        density[Shift_koor(idx_i      ,idx_j_cl   ,idx_k_cl   ,Ngrid,0)] +=  particle_weight*d_x    *(1-d_y)*(1-d_z)  ;
        density[Shift_koor(idx_i_cl   ,idx_j_cl   ,idx_k_cl   ,Ngrid,0)] +=  particle_weight*(1-d_x)*(1-d_y)*(1-d_z)  ;
        one_par = d_x*d_y*d_z  + (1-d_x)*d_y*d_z +\
        d_x*(1-d_y)*d_z   + d_x*d_y*(1-d_z) + (1-d_x)*(1-d_y)*d_z  +\
        (1-d_x)*d_y*(1-d_z)  + d_x*(1-d_y)*(1-d_z)+(1-d_x)*(1-d_y)*(1-d_z);
        nu_par += one_par;
        if (fabs(nu_par-idx_f) >10 )  printf("both nu_par - idx_f %f,%d,%f,%f,%f,%f,%f\n",nu_par,idx_f,x_p,y_p,z_p,particle_weight,one_par);
        if((d_x<0)||(d_x>1)) printf("something went wrong CIC %f\t%f\t%f\t%f\n",x_p,x_c,d_x,Box/2);
        if((d_y<0)||(d_y>1)) printf("something went wrong CIC %f\t%f\t%f\t%f\n",y_p,y_c,d_y,Box/2);
        if((d_z<0)||(d_z>1)) printf("something went wrong CIC %f\t%f\t%f\t%f\n",z_p,z_c,d_z,Box/2);
      }

      printf("...Done. \n Correct CIC Number of particles asigned  %f and actual number of particles  %d in box of size %f\n",nu_par,n_gal,Box);
      return density;
}
double Weight_CIC(float chi1,float chi2,struct magni *mag_arr,struct survey_functions *surv_funcs,int mode){
      /*COmpute weight of particle for liger densities*/
      /*Upadated for the new equation Next is CIC*/
      double b_z_bar,n_g_bar_red,weight,q_mag,evol_b_red,z_bar,z_red, chi_particle, chi_particle_red,epsiolon;
      int idx_surv,idx_surv_red;
      if (mode <4){
        printf("Mode<4 Impossible \n");
        exit(0);
      }else{
        chi_particle = chi2;
        chi_particle_red = chi1;
        idx_surv_red = Find_cell_value(chi_particle_red,surv_funcs);
        // epsiolon = get_interpolated_surv_value(idx_surv_red,chi_particle_red,surv_funcs,3);
        z_bar = gsl_spline_eval(surv_splines.spline_red, chi_particle, surv_splines.acc_red);
        z_red = gsl_spline_eval(surv_splines.spline_red, chi_particle_red, surv_splines.acc_red);
        if(mode==4)       weight = mag_arr[0].mag_tot;
        else if(mode==5)  weight = mag_arr[0].mag_tot_obs;
        else if(mode==6)  weight = mag_arr[0].mag_tot;
        else if(mode==10) weight = mag_arr[0].mag_tot_obs;
        else if(mode==11) weight = mag_arr[0].mag_tot;
        else if(mode==12) weight = mag_arr[0].mag_tot_obs;
        else weight = (z_red-z_bar)/(z_red+1.);
        // else weight = epsiolon*(z_red-z_bar)/(z_red+1);
        }
      return weight;
}
int Find_cell_value(float chi,struct survey_functions *srv_funcs){
      /*Float binary search find index of element closest to chi in survey function strucutres
      needs to be modified for generic arrays */
      float chi_middle;
      int first=0,last=0,middle=0,idx_of_chi=0 ,old_mid=0;
      last=surv_number_pts-1;
      if (chi>srv_funcs[last].chi) return surv_number_pts-1;
      middle = (first+last)/2;
      chi_middle =srv_funcs[middle].chi;
      while (old_mid!=middle) {
        old_mid = middle;
        if (chi>chi_middle) {
          first = middle;
          middle = (first+last)/2;
          chi_middle =srv_funcs[middle].chi;
        }else if (chi<chi_middle){
          last = middle;
          middle = (first+last)/2;
          chi_middle =srv_funcs[middle].chi;
        }
      }
      idx_of_chi= old_mid;
      if (idx_of_chi>0) {while (srv_funcs[idx_of_chi].chi>chi) idx_of_chi--;}
      if ((srv_funcs[idx_of_chi+1].chi<chi)){while (srv_funcs[idx_of_chi+1].chi<chi) idx_of_chi++;}
      return idx_of_chi;
}
float *Create_galaxy_density(float *dens_real_w,float *dens_shift_w,float *dens_evol,float *dens_mag,struct survey_functions *surv_func,int Ngrid, float gridsize,float max_particle_dist,int mode,float q_value){
      int Ngrid_loop = Ngrid/2,Ngrid3d = Ngrid*Ngrid*Ngrid,idx_dens,idx_surv,zero_dens = 0,neg_dens = 0,nan_dens = 0;
      double n_g_bar,x_p,y_p,z_p,chi,grid_volume = gridsize*gridsize*gridsize,non_zero_ratio,b_lin,w_g,epsiolon=0,q_mag=0;
      float *gal_density;
      if (!(gal_density = (float *) calloc( Ngrid3d , sizeof(float)) )) {printf("Failed to allocate memory for the grid.\n"); exit(0);}
      for (int idx_x = 0; idx_x < Ngrid/2; idx_x++){
        x_p = (float) idx_x *gridsize;
        if (x_p > max_particle_dist+gridsize) { Ngrid_loop = idx_x + 3;
          break;
        }
      }
      if (Ngrid_loop>Ngrid/2){printf("How is that possible %d,%d\n",Ngrid,Ngrid_loop );exit(0);}
      // max_particle_dist -= 2*gridsize;

      printf("Creating galaxy density up to %f with %d^3 iterations and %f...",max_particle_dist,Ngrid_loop,q_value);fflush(stdout);
      for (int idx_x = -1*Ngrid_loop; idx_x < Ngrid_loop; idx_x++) {
        for (int idx_y = -1*Ngrid_loop; idx_y < Ngrid_loop; idx_y++){
          for (int idx_z = -1*Ngrid_loop; idx_z < Ngrid_loop; idx_z++){
            x_p = (float) (idx_x + 0.5) *gridsize;
            y_p = (float) (idx_y + 0.5) *gridsize;
            z_p = (float) (idx_z + 0.5) *gridsize;
            chi = sqrt(x_p*x_p + y_p*y_p + z_p*z_p);if (chi > max_particle_dist) {zero_dens++; continue;}
            idx_dens = Shift_koor(idx_x      ,idx_y      ,idx_z      ,Ngrid,1);
            idx_surv = Find_cell_value(chi,surv_func);
            n_g_bar = get_interpolated_surv_value(idx_surv,chi,surv_func,0);
            b_lin = get_interpolated_surv_value(idx_surv,chi,surv_func,1);
            w_g=n_g_bar/n_sim_bar;
            if (mode == 0) gal_density[idx_dens] = b_lin*w_g*dens_real_w[idx_dens]  - n_g_bar*grid_volume*(b_lin-1);
            else{
              epsiolon = get_interpolated_surv_value(idx_surv,chi,surv_func,3);
              if (q_value <0 ) q_mag = get_interpolated_surv_value(idx_surv,chi,surv_func,2);
              else q_mag = q_value;

              if (dens_shift_w[idx_dens]>0.0){
                gal_density[idx_dens] = (b_lin - 1)*(w_g*dens_real_w[idx_dens] - n_g_bar*grid_volume) + w_g*dens_shift_w[idx_dens] \
                                        + n_g_bar*grid_volume*q_mag*dens_mag[idx_dens]/dens_shift_w[idx_dens] - n_g_bar*grid_volume*q_mag \
                                        + n_g_bar*grid_volume*epsiolon*dens_evol[idx_dens]/dens_shift_w[idx_dens];
                }
              else gal_density[idx_dens] = (b_lin - 1)*(w_g*dens_real_w[idx_dens] - n_g_bar*grid_volume) + w_g*dens_shift_w[idx_dens];
            }
            if ((chi<chi_neg)||(chi > max_particle_dist)||(gal_density[idx_dens]<0)) gal_density[idx_dens]=0.0;
            if ((gal_density[idx_dens]<0)&&(chi<max_particle_dist)) {neg_dens++;gal_density[idx_dens]=0;}else if (isnan(gal_density[idx_dens])) nan_dens++; else if (gal_density[idx_dens] == 0) zero_dens++;

          }
        }
      }
      non_zero_ratio = ((float)(Ngrid3d-neg_dens - nan_dens - zero_dens))/Ngrid3d;
      printf("....Galaxy Density formed with %d Nan, %d negative, %d zero elements. %f of the elements have values.\n",nan_dens, neg_dens,zero_dens,non_zero_ratio );fflush(stdout);
      fflush(stdout);
      return gal_density;
}
int Shift_koor(int idx_x,int idx_y,int idx_z,int Ngrid,int mode){
      int dum;
      if(mode){
        idx_x += Ngrid/2;
        idx_y += Ngrid/2;
        idx_z += Ngrid/2;
      }
      if (idx_x >= Ngrid){dum = idx_x/Ngrid; idx_x = idx_x - Ngrid*dum;}
      if (idx_y >= Ngrid){dum = idx_y/Ngrid; idx_y = idx_y - Ngrid*dum;}
      if (idx_z >= Ngrid){dum = idx_z/Ngrid; idx_z = idx_z - Ngrid*dum;}
      if ((idx_x>=Ngrid)||(idx_y>=Ngrid)||(idx_z>=Ngrid)||(idx_x<0)||(idx_y<0)||(idx_z<0)) {printf("ERROR:Shift_koor:%i (idx_x=%i/idx_y=%i/idx_z=%i/MAX=%i)\n",idx_z+Ngrid*(idx_y+Ngrid*idx_x),idx_x,idx_y,idx_z,Ngrid*Ngrid*Ngrid); exit(0);}
      return idx_z + Ngrid*(idx_y + Ngrid*idx_x);
    }
double get_interpolated_surv_value (int idx_surv,float chi, struct survey_functions *srv_funcs,int n_b_q){
      /*interpolates survey function to chi
      n_b_q: 0:n_z,1:b_lin,2:q_mag*/
      double result;
      if (USESPLINE && (chi > 50)){
        if(n_b_q == 0){
          result = gsl_spline_eval(surv_splines.spline_n_z, chi, surv_splines.acc_n_z);
        }else if (n_b_q == 1){
          result = gsl_spline_eval(surv_splines.spline_b_lin, chi, surv_splines.acc_b_lin);
        }else if (n_b_q == 2){
          result = gsl_spline_eval(surv_splines.spline_q_mag, chi, surv_splines.acc_q_mag);
        }else if (n_b_q == 3){
          result = gsl_spline_eval(surv_splines.spline_e_bias, chi, surv_splines.acc_e_bias);
        }
        return result;
      }else{//The else is not needed, I just think it's more readable
      if(n_b_q == 0){
        if (idx_surv == 0)                     result = Cubic_spline(srv_funcs[idx_surv+2].n_z,srv_funcs[idx_surv+1].n_z,srv_funcs[idx_surv].n_z,srv_funcs[idx_surv].n_z,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv].chi,chi);
        else if(idx_surv == surv_number_pts-2) result = Cubic_spline(srv_funcs[idx_surv+1].n_z,srv_funcs[idx_surv+1].n_z,srv_funcs[idx_surv].n_z,srv_funcs[idx_surv-1].n_z,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
        else if(idx_surv == surv_number_pts-1) result = srv_funcs[idx_surv].n_z;
        else                                   result = Cubic_spline(srv_funcs[idx_surv+2].n_z,srv_funcs[idx_surv+1].n_z,srv_funcs[idx_surv].n_z,srv_funcs[idx_surv-1].n_z,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
      }else if (n_b_q == 1){
        if (idx_surv == 0)                     result = Cubic_spline(srv_funcs[idx_surv+2].b_lin,srv_funcs[idx_surv+1].b_lin,srv_funcs[idx_surv].b_lin,srv_funcs[idx_surv].b_lin,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv].chi,chi);
        else if(idx_surv == surv_number_pts-2) result = Cubic_spline(srv_funcs[idx_surv+1].b_lin,srv_funcs[idx_surv+1].b_lin,srv_funcs[idx_surv].b_lin,srv_funcs[idx_surv-1].b_lin,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
        else if(idx_surv == surv_number_pts-1) result = srv_funcs[idx_surv].b_lin;
        else                                   result = Cubic_spline(srv_funcs[idx_surv+2].b_lin,srv_funcs[idx_surv+1].b_lin,srv_funcs[idx_surv].b_lin,srv_funcs[idx_surv-1].b_lin,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
      }else if (n_b_q == 2){
        if (idx_surv == 0)                     result = Cubic_spline(srv_funcs[idx_surv+2].q_mag,srv_funcs[idx_surv+1].q_mag,srv_funcs[idx_surv].q_mag,srv_funcs[idx_surv].q_mag,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv].chi,chi);
        else if(idx_surv == surv_number_pts-2) result = Cubic_spline(srv_funcs[idx_surv+1].q_mag,srv_funcs[idx_surv+1].q_mag,srv_funcs[idx_surv].q_mag,srv_funcs[idx_surv-1].q_mag,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
        else if(idx_surv == surv_number_pts-1) result = srv_funcs[idx_surv].q_mag;
        else                                   result = Cubic_spline(srv_funcs[idx_surv+2].q_mag,srv_funcs[idx_surv+1].q_mag,srv_funcs[idx_surv].q_mag,srv_funcs[idx_surv-1].q_mag,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
      }else if (n_b_q == 3){
        if (idx_surv == 0)                     result = Cubic_spline(srv_funcs[idx_surv+2].e_evl,srv_funcs[idx_surv+1].e_evl,srv_funcs[idx_surv].e_evl,srv_funcs[idx_surv].e_evl,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv].chi,chi);
        else if(idx_surv == surv_number_pts-2) result = Cubic_spline(srv_funcs[idx_surv+1].e_evl,srv_funcs[idx_surv+1].e_evl,srv_funcs[idx_surv].e_evl,srv_funcs[idx_surv-1].e_evl,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
        else if(idx_surv == surv_number_pts-1) result = srv_funcs[idx_surv].e_evl;
        else                                   result = Cubic_spline(srv_funcs[idx_surv+2].e_evl,srv_funcs[idx_surv+1].e_evl,srv_funcs[idx_surv].e_evl,srv_funcs[idx_surv-1].e_evl,srv_funcs[idx_surv+2].chi,srv_funcs[idx_surv+1].chi,srv_funcs[idx_surv].chi,srv_funcs[idx_surv-1].chi,chi);
      }
      return result;
    }
}
double Cubic_spline(double fieldfarer,double fieldfar,double fieldclose,double fieldcloser,double distfarer,double distfar,double distclose,double distcloser,double dist){
    double pg,g,p,c,d,dfdxc,dfdxf;

    g=distfar-distclose;
    p=(dist-distclose);
    if ((p < 0)||(p > g)) {
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
void Radial_bins(float *galaxy_density,struct survey_functions *srv_funcs,float z_min, float z_max,float grid_size,int Ngrid){
    double x_pos,y_pos,z_pos,chi,chi_min,chi_max,av=0,av_gal=0.0,grid_volume=grid_size*grid_size*grid_size,cell_nu=0,min_dens=2.0,max_dens=-1;
    int idx_dens,Ngrid_loop=Ngrid/2,idx_surv_gal;
    chi_min = calc_comoving_distance(z_min);
    chi_max = calc_comoving_distance(z_max);
    printf("CONSISTENCY CHECK...." );fflush(stdout);
    for (int idx_x = -1*Ngrid_loop; idx_x < Ngrid_loop; idx_x++) {
      for (int idx_y = -1*Ngrid_loop; idx_y < Ngrid_loop; idx_y++){
        for (int idx_z = -1*Ngrid_loop; idx_z < Ngrid_loop; idx_z++){
          x_pos = (float) (float) (idx_x + 0.5) *grid_size;
          y_pos = (float) (float) (idx_y + 0.5) *grid_size;
          z_pos = (float) (float) (idx_z + 0.5) *grid_size;
          chi = sqrt(x_pos*x_pos + y_pos*y_pos + z_pos*z_pos);
          idx_dens=Shift_koor(idx_x,idx_y,idx_z,Ngrid,1);
          if((chi<chi_max)&&(chi>chi_min)){
            // if(galaxy_density[idx_dens]<0) printf("negativ %f,%f,(%d,%d,%d) \n",galaxy_density[idx_dens],chi,idx_x,idx_y,idx_z );
            idx_surv_gal=Find_cell_value(chi,srv_funcs);
            av_gal+=srv_funcs[idx_surv_gal].n_z;//*(srv_funcs[idx_surv_gal].b_lin-1);
            if(galaxy_density[idx_dens]>max_dens) max_dens = galaxy_density[idx_dens];
            if((galaxy_density[idx_dens]<min_dens)&&(galaxy_density[idx_dens]>0)) min_dens = galaxy_density[idx_dens];
            av+=galaxy_density[idx_dens];
            cell_nu+=1.0;
          }
        }
      }
    }
    av/=cell_nu;
    av_gal/=cell_nu;

    printf("Done.\n Angular test (Average input density, n_g, ratio) in range (z_i,z_f) ,(chi_min,chi_max),(min ,max)\n %f,%f,%f,(%f,%f),(%f,%f),(%.8f,%.8f)\n",av,av_gal*grid_volume,av/(av_gal*grid_volume) ,z_min,z_max,chi_min,chi_max,min_dens, max_dens);
    if (fabs(av/(av_gal*grid_volume)-1)>2e-2) printf("WARNING:Ratio is too far from one\n");
    printf("___________________________________________________________");
    fflush(stdout);
}
void ADD_shot(float *corr_dens,int Ngrid,float box_length){
      int Ngrid3d = Ngrid*Ngrid*Ngrid;
      float dum;
      float const_norm = 1.;
      const gsl_rng_type * T;
      gsl_rng * r;
      unsigned int N_possion;
      gsl_rng_env_setup();
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);
      printf("Possion sampling each cell\n");
      fflush(stdout);
      for (int idx_i = 0; idx_i < Ngrid3d; idx_i++)
        {
          dum = corr_dens[idx_i];
          N_possion= gsl_ran_poisson (r, dum*const_norm);
          corr_dens[idx_i]= ((float) N_possion)/const_norm;
        }
      printf("\tDone\n");
      fflush(stdout);
      gsl_rng_free (r);
    }
