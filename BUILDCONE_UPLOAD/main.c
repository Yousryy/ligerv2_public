#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "src/vardef.h"
#include "src/IO.h"
#include "src/create_galaxy_density.h"
#include "src/Input_read.h"
int main(int argc, char const *argv[]) {
    char param_file[200];
    sprintf(param_file, argv[1]);
    int params;
    fflush(stdout);
    params=read_ParaFile(param_file);
    Read_simulation_parameters(run.sim_path);
    Build_Cone(run.liger_path,run.sim_path,run.surv_func_path,run.save_dir_path,run.Ngrid);
    // if (run.COMPUTE_SPECTRA = 1)

    fflush(stdout);
  return 0;
}
