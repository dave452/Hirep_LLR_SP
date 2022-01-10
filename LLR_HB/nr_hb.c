/****************************************************************************
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"
#include "spectrum.h"
#include "llr_hb_utils.h"
#include "cinfo.c"
//#include "wilsonflow.h"
#include "setup.h"

/* LLR parameters */
typedef struct _input_llr {
  char make[256];
  int nmc,nth,it;
  double starta,S0,dS;
  /* for the reading function */
  input_record_t read[8];
} input_llr;


#define init_input_llr(varname) \
  { \
  .read={\
    {"make llr iterations", "llr:make = %s", STRING_T, &((varname).make)}, \
    {"Number of MC steps per RM iteration ", "llr:nmc = %d", INT_T, &((varname).nmc)}, \
    {"Number of MC therm steps per RM iteration", "llr:nth = %d", INT_T, &((varname).nth)}, \
    {"Initial a", "llr:starta = %lf", DOUBLE_T, &((varname).starta)}, \
    {"Robbins Monro startint iteration", "llr:it = %d", INT_T, &((varname).it)}, \
    {"Cental action", "llr:S0 = %lf", DOUBLE_T, &((varname).S0)}, \
    {"Delta S", "llr:dS = %lf", DOUBLE_T, &((varname).dS)}, \
    {NULL, NULL, 0, NULL}				\
    }\
}
 
input_llr llr_var=init_input_llr(llr_var);

pg_flow flow=init_pg_flow(flow);

char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";

static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: [-i <input file>] [-o <output file>] [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}



int main(int argc,char *argv[]) {
  struct timeval startmain, endmain, etimemain; /* //for trajectory timing */
  gettimeofday(&startmain,0);
  char sbuf[128];
  int i;
  read_cmdline(argc,argv);
  
  /* setup process communications */
  setup_process(&argc,&argv);
  
  /* read global variables file */
  read_input(glb_var.read,input_filename);
  
  setup_replicas();
  
  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    FILE* stderrp;
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    stderrp=freopen(error_filename,"w",stderr);
    error(stderrp==NULL,1,"main [hmc.c]",
	  "Cannot redirect the stderr");
  }
  
  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS);
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);
  //lprintf("MAIN",0,"SVN Revision: %d\n", CI_svnrevision);

  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));
  
  /* setup lattice geometry */
  //if (geometry_init() == 1) { finalize_process(); lprintf("TEST",0,"Geometry_init: \n"); return 0; }
  //geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */ 
  //lprintf("TEST",0,": geometry_mpi_eo \n");
  /* setup random numbers */
  read_input(rlx_var.read,input_filename);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

  if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0')
  {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
    read_ranlxd_state(rlx_var.rlxd_state);
  }

#ifdef GAUGE_SUN
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Default gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* read input for llr update */
  read_input(llr_var.read,input_filename);

  lprintf("MAIN",0,"LLR nunber of mc steps per RM: %d\n",llr_var.nmc);
  lprintf("MAIN",0,"LLR nunber of therm steps per RM %d\n",llr_var.nth);
  lprintf("MAIN",0,"LLR Initial a %f\n",llr_var.starta);
  lprintf("MAIN",0,"LLR RM start value iteration %d\n",llr_var.it);
  lprintf("MAIN",0,"LLR S0 Central action %f\n",llr_var.S0);
  lprintf("MAIN",0,"LLR Delta S %f\n",llr_var.dS);
 
  /* Init Monte Carlo */

  init_mc(&flow, input_filename);
  //lprintf("MAIN",0,"MVM during HMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

  //for(int j=0;j<100;++j) {
  //       update(llr_var.starta,1,0);
  //       //update(2.4,1,0);
  //       //printf(" E= %1.8e, %f\n", avr_plaquette(), avr_plaquette()*6.0*GLB_VOLUME);
  //}
  
  //double E = avr_plaquette()*GLB_VOLUME*6.;
  init_robbinsmonro(llr_var.nmc,llr_var.nth,llr_var.starta,llr_var.it,llr_var.dS,llr_var.S0);
  

  for(int j=0;j<flow.rmrestart;++j) {
    
    restart_robbinsmonro();  
    
    for (i=0;i<flow.therm;++i){
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"Starting thermalization, thermalization steps = %d\n",flow.therm);
      lprintf("MAIN",0,"Thermalization #%d...\n",i);
      gettimeofday(&start,0);
      
      thermrobbinsmonro();
      
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Thermalization sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);    
    }
    
    //lprintf("MAIN",0,"Thermalization done.\n");
    

    for(i=flow.start;i<flow.end;++i) {
      
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"Trajectory #%d...\n",i);
      gettimeofday(&start,0);
      
     newtonraphson();
    //  
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Robbins Monro sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      lprintf("MAIN",0,"Plaq a fixed %lf \n",avr_plaquette());    
      lprintf("MAIN",0,"<a_rho(%d,%d,%lf)>= %f\n",j,i,getS0(),get_llr_a());  
    }
    lprintf("MAIN",0,"Robins Monro update done.\n");
  /* save final configuration */
  save_conf(&flow, (flow.end-flow.start)*flow.rmrestart);
  ///* Only save state if we have a file to save to */
  if(rlx_var.rlxd_state[0]!='\0') {
    lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
     write_ranlxd_state(rlx_var.rlxd_state);
  }
  gettimeofday(&endmain,0);
  timeval_subtract(&etimemain,&endmain,&startmain);
  
  lprintf("MAIN",0,"Total simulation time =[%ld sec %ld usec]\n",etimemain.tv_sec,etimemain.tv_usec);
  /* finalize Monte Carlo */
  end_mc();
  
  /* close communications */
  finalize_process();
  
  return 0;
  
}
