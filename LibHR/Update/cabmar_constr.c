/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File cabmar.c
*
* Cabibbo-Marinari rotations
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "update.h"

#include "observables.h"
#include "global.h"

#ifndef WITH_QUATERNIONS

static inline void rotate(suNg_vector *pu1, suNg_vector *pu2, double *s)
{
  int i;
  complex z1,z2;
  complex *cu1, *cu2;
  
  cu1 = &((*pu1).c[0]);
  cu2 = &((*pu2).c[0]);
  
  for (i=0; i<NG; ++i) {
    z1=s[0]*creal(*cu1)-s[1]*cimag(*cu2)+s[2]*creal(*cu2)-s[3]*cimag(*cu1)+I*(s[0]*cimag(*cu1)+s[1]*creal(*cu2)+s[2]*cimag(*cu2)+s[3]*creal(*cu1));
    z2=s[0]*creal(*cu2)-s[1]*cimag(*cu1)-s[2]*creal(*cu1)+s[3]*cimag(*cu2)+I*(s[0]*cimag(*cu2)+s[1]*creal(*cu1)-s[2]*cimag(*cu1)-s[3]*creal(*cu2));
    (*cu1) = z1;
    (*cu2) = z2;
    ++cu1;
    ++cu2;
  }
}



static inline void wmatrix(suNg_vector *pu1, suNg_vector *pu2, suNg_vector *pv1, suNg_vector *pv2, double *w)
{
	double prod1,prod2;
  _vector_prod_re_g(prod1,*pu1,*pv1);
	_vector_prod_re_g(prod2,*pu2,*pv2);
	w[0] = prod1+prod2;
  _vector_prod_im_g(prod1,*pu1,*pv2);
	_vector_prod_im_g(prod2,*pu2,*pv1);
	w[1] = prod1+prod2;
  _vector_prod_re_g(prod1,*pu2,*pv1);
	_vector_prod_re_g(prod2,*pu1,*pv2);
	w[2] = prod1-prod2;
  _vector_prod_im_g(prod1,*pu1,*pv1);
	_vector_prod_im_g(prod2,*pu2,*pv2);
	w[3] = prod1-prod2;
}

#endif

#ifndef GAUGE_SUN
#ifdef LLRHB
#error LLR HB only available for SUN
#endif
#else

void cabmar_constrained(double beta,suNg *u,suNg *v,int type, double * E, double Emin, double Emax)
{
	double xmin, xmax, Eold, Enew, sold, snew;
	suNg W1;
#ifdef WITH_QUATERNIONS
  double wsq,rho,fact;
	suNg s,w,r;
  
	/*w=u*v^+ il dagger c'è perché la staple
	 * in Hirep è definita percorsa al contrario*/
	_suNg_times_suNg_dagger(w,*u,*v);
	
	Eold = w.c[0];
	
	/*wsq=w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];*/
	_suNg_sqnorm(wsq,w); //sqnorm gives twice what we need
  wsq*=0.5;
  
	if ((beta*beta*wsq)>1.0e-28) {
		if (type==1) {
			_suNg_trace_re(fact,w);
			fact/=wsq;
			_suNg_mul(s,fact,*v);
			/* _suNg_sub(*u,s,*u); */
      _suNg_minus(*u,*u);
      _suNg_add_assign(*u,s);
		} else {
			fact=sqrt(wsq);
			rho=beta*fact;
			xmax = (-*E + Emax + Eold)/fact;
			xmin = (-*E + Emin + Eold)/fact;
			xmin=(xmin>-1.)?xmin:-1.;
			xmax=(xmax<1.)?xmax:1.;
			//printf("========================\n");
			//printf("Emeas = %f, %f\n", avr_plaquette()*6.0*GLB_VOLUME, avr_plaquette());
			//printf("E = %f, Eold= %f, Emin= %f, Emax= %f\n", *E, Eold, Emin, Emax);
			//printf("xmin = %f, xmax = %f\n", xmin, xmax);
			random_su2_constrained(rho,&r.c[0],xmin,xmax);
			//random_su2(rho,&r.c[0]);
      
			fact=1./fact;

      			_suNg_times_suNg(*u,r,*v);
			_suNg_mul(*u,fact,*u);
			_suNg_times_suNg_dagger(w,*u,*v);
			//Snew = 1.-w.c[0];
			Enew = w.c[0];
			*E += (Enew-Eold);
		}
	} else ;//random_su2_constrained(0.0,&(u->c[0]),-1.,1.);
	
  
#else
  int i,j;
  double b,bsq,wsq,rho,fact,r[4],w[4],s[4];
  const double invng = 1. / (double) NG;
  
  
  suNg_vector *pu1=(suNg_vector*)(u);
  suNg_vector *pv1=(suNg_vector*)(v);
  
  b=invng*beta;
  bsq=b*b;
  
  for (i=0; i<NG; ++i) {
    suNg_vector *pu2 = pu1 + 1;
    suNg_vector *pv2 = pv1 + 1;
    for (j=i+1; j<NG; ++j) {
      

      wmatrix(pu1, pu2, pv1, pv2, w);
      wsq=w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];
      Eold = invng*w[0]; 
      fact=sqrt(wsq);
      rho=b*fact;
      //giuste!
      xmax = NG*(-*E+Emax+Eold)/(fact);
      xmin = NG*(-*E+Emin+Eold)/(fact);
      //
      //xmax = -NG*(-*S+Smin+sold)/(fact);
      //xmin = -NG*(-*S+Smax+sold)/(fact);
      
      //printf("xmin = %f, xmax = %f\n", xmin, xmax);
      xmin=(xmin>-1.)?xmin:-1.;
      xmax=(xmax<1.)?xmax:1.;
      //printf("==============================\n");
      //printf("Emeas = %f, %f\n", avr_plaquette()*6.0*GLB_VOLUME, avr_plaquette());
      //printf("E = %f, Eold = %f, Emin= %f, Emax= %f\n", *E, Eold, Emin, Emax);
      //printf("xmin = %f, xmax = %f\n", xmin, xmax);
      
      if ((bsq*wsq)>1.0e-28 ) {
        if (type==1) {
          fact=(w[0]+w[0])/wsq;
          s[0]=fact*w[0]-1.0;
          s[1]=fact*w[1];
          s[2]=fact*w[2];
          s[3]=fact*w[3];
        } else {
          //Sold = -invng*w[0]; 
          //printf("Emeas = %f, %f\n", avr_plaquette()*6.0*GLB_VOLUME, avr_plaquette());
          
          random_su2_constrained(rho,r, xmin,xmax);
          
          fact=1.0/fact;
          s[0]=fact*(r[0]*w[0]-r[1]*w[1]-r[2]*w[2]-r[3]*w[3]);
          s[1]=fact*(r[1]*w[0]+r[0]*w[1]-r[2]*w[3]+r[3]*w[2]);
          s[2]=fact*(r[2]*w[0]+r[0]*w[2]-r[3]*w[1]+r[1]*w[3]);
          s[3]=fact*(r[3]*w[0]+r[0]*w[3]-r[1]*w[2]+r[2]*w[1]);
	  rotate(pu1,pu2,s);
        }
      } else ;//random_su2(0.0,s);
      
      //printf("Emeas = %f, %f\n", avr_plaquette()*6.0*GLB_VOLUME, avr_plaquette());
      wmatrix(pu1, pu2, pv1, pv2, w);
      //_suNg_times_suNg_dagger(W1,*u,*v);
      //_suNg_trace_re(Snew,W1); 
      //Snew *= invng;
      Enew = invng*w[0];

      //printf("E = %f, Enew- Eold = %f, Emin= %f, Emax= %f\n", *S, Snew-Sold, Smin, Smax);
      *E += (Enew - Eold); 
      //printf("Emeas = %f, %f\n", avr_plaquette()*6.0*GLB_VOLUME, avr_plaquette());
      ++pu2; ++pv2;
    }
    ++pu1; ++pv1;
  }
#endif
}
#endif
