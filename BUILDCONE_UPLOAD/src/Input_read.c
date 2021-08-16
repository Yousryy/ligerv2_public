#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./Input_read.h"
#include "./vardef.h"


#define NUMPARAS 20
#define TEXTPARA 0
#define FLOATPARA 1
#define INTPARA 2
#define TRUTHPARA 3
#define INTARRAYPARA 4
struct Global_VAR run;
int cmpfunc (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}
int read_ParaFile(char *fname){
  char Paras[NUMPARAS][100],ReadLine[500],ReadParaBuf[200],RV[11][200];
  FILE *Pfile;
  int ParaType[NUMPARAS],Npara,i,Nread;
  void *ParaAdr[NUMPARAS];
  Npara=0;
  sprintf(Paras[Npara],"Light_cone_path");
  ParaType[Npara]=TEXTPARA;
  ParaAdr[Npara]=run.liger_path;
  ++Npara;

  sprintf(Paras[Npara],"simulation_path");
  ParaType[Npara]=TEXTPARA;
  ParaAdr[Npara]=run.sim_path;
  ++Npara;

  sprintf(Paras[Npara],"Survey_func_path");
  ParaType[Npara]=TEXTPARA;
  ParaAdr[Npara]=run.surv_func_path;
  ++Npara;

  sprintf(Paras[Npara],"Save_base_dir");
  ParaType[Npara]=TEXTPARA;
  ParaAdr[Npara]=run.save_dir_path;
  ++Npara;

  int *wanted_dens;
  sprintf(Paras[Npara],"Box_size_factor");
  ParaType[Npara]=FLOATPARA;
  ParaAdr[Npara]=&run.factor_box_size;
  ++Npara;

  sprintf(Paras[Npara],"Q_const");
  ParaType[Npara]=FLOATPARA;
  ParaAdr[Npara]=&run.Q_const;
  ++Npara;

  sprintf(Paras[Npara],"NGRID");
  ParaType[Npara]=INTPARA;
  ParaAdr[Npara]=&run.Ngrid;
  ++Npara;

  sprintf(Paras[Npara],"wanted_dens");
  ParaType[Npara]=INTARRAYPARA;
  ParaAdr[Npara]=&run.dens_to_compute[0];
  ++Npara;

  sprintf(Paras[Npara],"SN");
  ParaType[Npara]=INTPARA;
  ParaAdr[Npara]=&run.SN;
  ++Npara;
  if((Pfile=fopen(fname, "r"))) {
    while (fgets(ReadLine, sizeof(ReadLine),Pfile) != NULL){
      Nread=sscanf(ReadLine,"%s%s%s%s%s%s%s%s%s%s%s%s%s",ReadParaBuf,RV[0],RV[1],RV[2],RV[3],RV[4],RV[5],RV[6],RV[7],RV[8],RV[9],RV[10],RV[11]);
      if (Nread < 2) continue;
      if (ReadParaBuf[0]==*"%") continue;
      if (ReadParaBuf[0]==*"#") continue;
      for (i=0;i<Npara;++i) if (strcmp(ReadParaBuf, Paras[i]) == 0) break;
      if (i==Npara) {
        printf("	Ignoring unkown parameter %s\n",ReadParaBuf);
        continue;
      }
      switch (ParaType[i]) {
        case FLOATPARA:
        ParaType[i]=-1;
        if (sscanf(RV[0],"%f",((float *)ParaAdr[i]))!=1){
          printf("	Error reading %s, should be an floating-point value!\n",ReadParaBuf);
          exit(0);
        }
        // printf("%s\n",RV[0] );
        break;

        case INTARRAYPARA:
        if (Nread > 11) {
          printf("	Expect less numerical values for %s\n",ReadParaBuf);
          exit(0);
        }
        run.nu_of_cones = Nread-1;
        if (!(run.dens_to_compute=malloc(run.nu_of_cones*sizeof(int)))) {printf("\nAllocation of memory failed.\n"); exit(0);}
        ParaType[i]=-1;
        for (int idx= 0; idx < run.nu_of_cones; idx++) {
          if (sscanf(RV[idx],"%i",&run.dens_to_compute[idx])!=1){
            printf("	Error reading %s, should be an int-point value!\n",ReadParaBuf);
            exit(0);
          }
        }
        qsort(&run.dens_to_compute[0], run.nu_of_cones, sizeof(int), cmpfunc);
        break;

        case TEXTPARA:
        ParaType[i]=-1;
        strcpy(ParaAdr[i], RV[0]); break;

        case INTPARA:
        ParaType[i]=-1;
        if (sscanf(RV[0],"%i",((int *)ParaAdr[i]))!=1){
          printf("	Error reading %s, should be an integer value!\n",ReadParaBuf);
          exit(0);
        } break;

        case -1:
        printf("	Parameter %s appears more then once!\n",ReadParaBuf);
        exit(0); break;
      }
    }
    fclose(Pfile);
  } else {
    printf("	Cannot read file with parameters from: %s\n",fname);
    exit(0);
  }

  Nread=0;
  for (i=0;i<Npara;++i) {
    if (ParaType[i]!=-1) printf("	Missing parameter: %s\n",Paras[i]);
    else ++Nread;
  }
  if (Nread!=Npara) {
    printf("%i parameters missing.\n",Npara-Nread);
    exit(0);
  }

  for (int idx = 0; idx < run.nu_of_cones; idx++) {
    if ((run.dens_to_compute[idx]<0)||(run.dens_to_compute[idx]>9)){printf("unavailable option for wanted_dens,%i\n",run.dens_to_compute[idx]);exit(0);}
    for (int idx2 = 0; idx2 < run.nu_of_cones; idx2++){
      if ((run.dens_to_compute[idx]==run.dens_to_compute[idx2])&&(idx!=idx2)){printf("repeated wanted_dens choice,%i\n",run.dens_to_compute[idx]);exit(0);}
    }
    if (idx>0){if(run.dens_to_compute[idx]<run.dens_to_compute[idx-1]){printf("Please sort wanted_dens \n");exit(0);}}
  }
  printf("The input densities indices are  \n " );
  for (int  idx_d = 0; idx_d < run.nu_of_cones; idx_d++) {
    printf("%d\t",run.dens_to_compute[idx_d]);
  }
  printf("\n");
  fflush(stdout);
  return Nread;
}

#undef NUMPARAS
#undef TEXTPARA
#undef FLOATPARA
#undef INTPARA
#undef TRUTHPARA
#undef INTARRAYPARA
