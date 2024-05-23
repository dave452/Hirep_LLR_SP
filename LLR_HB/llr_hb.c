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
#include "setup.h"
//#include "wilsonflow.h"


/* LLR parameters */
typedef struct _input_llr {
  char make[256];
  int nmc,nth,it, nfxa, sfreq_fxa, nhb, nor, it_freq, N_nr;
  double starta,S0,dS, Smin, Smax, db;
  /* for the reading function */
  input_record_t read[17];
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
    {"Maximum S value for all replicas", "llr:Smax = %lf", DOUBLE_T, &((varname).Smax)}, \
    {"Minimum S value for all replicas", "llr:Smin = %lf", DOUBLE_T, &((varname).Smin)}, \
    {"Number of fixed a steps ", "llr:nfxa = %d", INT_T, &((varname).nfxa)}, \
    {"Swap frequency for fixed a interations ", "llr:sfreq_fxa = %d", INT_T, &((varname).sfreq_fxa)}, \
    {"Number of heatbath steps per MC step ", "nhb = %d", INT_T, &((varname).nhb)}, \
    {"Number of over-relaxation steps per MC step ", "nor = %d", INT_T, &((varname).nor)}, \
    {"Suppresion factor increment frequency ", "llr:it_freq = %d", INT_T, &((varname).it_freq)},  \
    {"Number of intial NR iterations  ", "llr:N_nr = %d", INT_T, &((varname).N_nr)},  \
    {"Change in beta for annealing ", "llr:db = %lf", DOUBLE_T, &((varname).db)},  \
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
  
  setup_gauge_fields();

  /* read input for llr update */
  read_input(llr_var.read,input_filename);

  lprintf("MAIN",0,"LLR number of mc steps per RM: %d\n",llr_var.nmc);
  lprintf("MAIN",0,"LLR number of therm steps per RM %d\n",llr_var.nth);
  lprintf("MAIN",0,"LLR Initial a %f\n",llr_var.starta);
  lprintf("MAIN",0,"LLR RM start value iteration %d and its iteration frequency %d \n",llr_var.it, llr_var.it_freq);
  lprintf("MAIN",0,"LLR S0 Central action %f\n",llr_var.S0);
  lprintf("MAIN",0,"LLR Delta S %f\n",llr_var.dS);
  lprintf("MAIN",0,"Number of heatbath steps %d, number of overrelaxation steps %d\n",llr_var.nhb,llr_var.nor);
  lprintf("MAIN",0,"Number of fixed a steps %d\n",llr_var.nfxa);
  lprintf("MAIN",0,"Swap frequency for fixed a steps %d\n",llr_var.sfreq_fxa);
  lprintf("MAIN",0,"LLR Smin minimum action for all replicas %f\n",llr_var.Smin);
  lprintf("MAIN",0,"LLR Smax maximum action for all replicas %f\n",llr_var.Smax);
  lprintf("MAIN",0,"LLR Delta S %f\n",llr_var.dS);
  lprintf("MAIN",0,"LLR dB %f\n",llr_var.db);
  lprintf("MAIN",0,"LLR number of intial NR iterations %d\n",llr_var.N_nr);
#ifdef LLRHBPARALLEL
  lprintf("MAIN",0,"Compiled with domain decomposition \n");
#else
  lprintf("MAIN",0,"Compiled without domain decomposition \n");
#endif
/* Init Monte Carlo */
  int initial_it;
  init_mc(&flow, input_filename);
  if(flow.start < llr_var.it)
  {
      flow.end = llr_var.it + (flow.end - flow.start);
      flow.start = llr_var.it;
  }
  else
  {
      llr_var.it = flow.start;

  }
  //initial_it = llr_var.it;
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

  init_robbinsmonro(llr_var.nmc,llr_var.nth,llr_var.starta,llr_var.it,llr_var.dS,llr_var.S0,llr_var.sfreq_fxa, llr_var.Smin, llr_var.Smax,llr_var.nhb,llr_var.nor, llr_var.it_freq, llr_var.db);


  for(int j=0;j<flow.rmrestart;++j) {

    restart_robbinsmonro();
    struct timeval start, end, etime; /* //for trajectory timing */
    gettimeofday(&start,0);
    for (i=0;i<flow.therm;++i){
      thermrobbinsmonro();
      if (flow.therm > 20)
        {
            if (i % (flow.therm / 5) == 0)
                lprintf("MAIN", 0, "%d", ((i * 100) / flow.therm));
            else if (i % (flow.therm / 20) == 0)
                lprintf("MAIN", 0, ".");
        }
    }
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Thermalization done in [%ld sec %ld usec].\n", etime.tv_sec, etime.tv_usec);
    lprintf("MAIN",0,"Measured energy %lf, stated energy \n",avr_plaquette()*6.* GLB_VOLUME,getS0() );
    lprintf("MAIN", 0, "Initial NR iterations, N_nr: %d %d \n", llr_var.N_nr);
    for(i=0;i<llr_var.N_nr;++i) {

      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"NR Trajectory #%d...\n",i);
      gettimeofday(&start,0);
      newtonraphson();
     //Timing and output data
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Newton Raphson sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      lprintf("MAIN",0,"NR Plaq a fixed %lf \n",avr_plaquette());
      lprintf("MAIN",0,"NR <a_rho(%d,%d,%.9f)>= %.9f\n",j,i,getS0(),get_llr_a());
    }
    //initial_it = llr_var.it;
    lprintf("MAIN",0,"Newton Raphson update done.\n");
    lprintf("MAIN", 0, "flow.start: %d, flow.end: %d, llr_var.it: %d \n", flow.start,flow.end, llr_var.it);
    for(i=flow.start;i<flow.end;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */
      lprintf("MAIN",0,"Trajectory #%d...\n",i);
      gettimeofday(&start,0);
      robbinsmonro();
     //Timing and output data
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Robbins Monro sequence #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      lprintf("MAIN",0,"Plaq a fixed %lf \n",avr_plaquette());
      lprintf("MAIN",0,"<a_rho(%d,%d,%.9f)>= %.9f\n",j,i,getS0(),get_llr_a());
    }
    lprintf("MAIN",0,"Robins Monro update done.\n");
    for(i=0;i<llr_var.nfxa;++i) {
      struct timeval start, end, etime; /* //for trajectory timing */

      llr_fixed_a_update();

      //lprintf("MAIN",0,"Obs measure for fixed E=%f dE=%f a=%f T\n",getS0(),getdS(),get_llr_a());
       }
  }
  /* save final configuration */
  save_conf(&flow, flow.end - 1);
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
