/*************************************************************************** \
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
 *
 * File update.c
 *
 * Update programs
 *
 *******************************************************************************/

#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"
#include "communications.h"
#include "random.h"
#include "observables.h"
#define PI 3.141592653589793238462643383279502884197
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static int *dyn_gauge=NULL;

//void project_gauge_field(void)
//{
//  _MASTER_FOR(&glattice,ix) {
//    project_to_suNg(pu_gauge(ix,0));
//    project_to_suNg(pu_gauge(ix,1));
//    project_to_suNg(pu_gauge(ix,2));
//    project_to_suNg(pu_gauge(ix,3));
//  }
//  
//  start_gf_sendrecv(u_gauge);
//} 

#if defined(BASIC_SF) || defined(ROTATED_SF)
static void g_up_Dirichlet_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BASIC_SF) || defined(ROTATED_SF) || defined(BC_T_MIXED)
static void g_dn_Dirichlet_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(1,ix,iy,iz);
	  dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN) || defined(BC_T_MIXED)
static void g_up_open_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(T-1,ix,iy,iz);
	  dyn_gauge[index*4]=0;
	}
  }
}
#endif

#if defined(BC_T_OPEN) 
static void g_dn_open_BCs() {
  int ix,iy,iz,index;
  
  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
	  index=ipt(0,ix,iy,iz);
	  dyn_gauge[index*4]=dyn_gauge[index*4+1]=dyn_gauge[index*4+2]=dyn_gauge[index*4+3]=0;
	}
  }
}
#endif

static void free_hb_boundary() {
  if (dyn_gauge!=NULL) {
    free(dyn_gauge);
    dyn_gauge = NULL;
  }
}

static void init_hb_boundary() {
  dyn_gauge = malloc(sizeof(*dyn_gauge)*glattice.gsize_gauge*4);
  atexit(&free_hb_boundary); //register cleanup function at exit
  
  for(int i=0;i<glattice.gsize_gauge*4;i++) dyn_gauge[i]=1;
#if defined(BASIC_SF) || defined(ROTATED_SF)
  g_up_Dirichlet_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_MIXED
  g_up_open_BCs();
  g_dn_Dirichlet_BCs();
#endif
#ifdef BC_T_OPEN
  g_up_open_BCs();
  g_dn_open_BCs();
#endif
}






static void update_all(double beta,int type, double * S, double Smin, double Smax)
{
  static int count=PROJECT_INTERVAL;
  
  if (count>=PROJECT_INTERVAL) {
    project_gauge_field();
    count=0;
  }
  ++count;
 //printf("boh= %d\n", glattice.local_master_pieces);
 //printf("boh= %d\n", glattice.master_start[0]);
 //printf("boh= %d\n", glattice.master_end[1]);
    suNg v;
//lprintf("MAIN",10,"Emin, Emax, E = %f, %f, %f\n", Smin, Smax, *S);
    for(int mu=0;mu<4;mu++){
      for(int i=0;i<glattice.local_master_pieces;i++) {
        for(int j=glattice.master_start[i];j<=glattice.master_end[i];j++){
          if(dyn_gauge[j*4+mu]!=0){
            staples(j,mu,&v);
	    //printf("E = %f, Emin= %f, Emax= %f, i=%d, j=%d\n", S,Smin,Smax,i,j);
            cabmar_constrained(beta,pu_gauge(j,mu),&v,type, S,Smin , Smax);
          }
    //printf("E = %f, Emin= %f, Emax= %f\n", *S, Smin, Smax);
        }
      }
    }
  
} 


void update_constrained(double beta,int nhb,int nor, double * S, double Smin, double Smax)
{
  if(dyn_gauge==NULL ) init_hb_boundary();
   
  //lprintf("MAIN", 10,"starting constrained update with a = %f\n", beta);
  for (int n=0;n<nhb;n++){
    //printf("update_all: E = %f, Emin= %f, Emax= %f\n", *S, Smin, Smax);
    update_all(beta,0, S, Smin, Smax);
  }

 // for (int n=0;n<nor;n++){
  //  update_all(beta,1);
 // }

  start_gf_sendrecv(u_gauge);
 
}

void anneal(double * S, double S0, double dS){
  suNg * u;
  suNg unew;
  //_suNg_unit(unew);
  suNg v,w,urand,w1;
  double k;
  if(dyn_gauge==NULL ) init_hb_boundary();
  for(;;){
  for(int mu=0;mu<4;mu++){
    for(int i=0;i<glattice.local_master_pieces;i++) {
      for(int j=glattice.master_start[i];j<=glattice.master_end[i] ;j++){
        if(dyn_gauge[j*4+mu]!=0){
            staples(j,mu,&v);
            random_suNg_epsilon(&urand,0.2);
            u = pu_gauge(j,mu);
            _suNg_times_suNg(unew,urand,*u);
            _suNg_sub_assign(unew,*u);
            _suNg_times_suNg_dagger(w,unew,v);
            _suNg_trace_re(k,w);
            k /= NG;
            if( ((*S > S0) && (k<0)) || ( (*S < S0) && (k>0)) ) {
                *S +=k;
                _suNg_times_suNg(w1,urand,*u);
                _suNg_mul(*u,1.,w1);
		//printf(" change made: S = %1.8e \n", *S);
  		//printf( "S = %1.8e\n", (1.-avr_plaquette())*6.0*GLB_VOLUME);
  		//printf( "E = %1.8e\n", (avr_plaquette())*6.0*GLB_VOLUME);
            } else ;
            if( abs(*S-S0) < 0.25*dS ) return;
         }
        }
      }
    }
  }
}

