/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#ifdef SGI

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "WeICME.h"
#include "sgi.h"

ITG *irowsgi=NULL;
double *ausgi=NULL;

void sgi_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
                ITG *neq, ITG *nzs, ITG token){

  char *oocpath="/yatmp/scr1",*env;
  ITG i,j,k,l,*pointers=NULL,method;
  long long ndim;
  double ops=0,ooclimit=2000.;

  printf(" Factoring the system of equations using the sgi solver\n\n");

  env=getenv("CCX_OOC_MEM");
  if(env) ooclimit=atoi(env);

  ndim=*neq+*nzs;

  NNEW(pointers,ITG,*neq+1);
  NNEW(irowsgi,ITG,ndim);
  NNEW(ausgi,double,ndim);

  k=ndim;
  l=*nzs;

  if(*sigma==0.){
    pointers[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowsgi[--k]=irow[--l]-1;
	ausgi[k]=au[l];
      }
      pointers[i]=--k;
      irowsgi[k]=i;
      ausgi[k]=ad[i];
    }
  }
  else{
    pointers[*neq]=ndim;
    for(i=*neq-1;i>=0;--i){
      for(j=0;j<icol[i];++j){
	irowsgi[--k]=irow[--l]-1;
	ausgi[k]=au[l]-*sigma*aub[l];
      }
      pointers[i]=--k;
      irowsgi[k]=i;
      ausgi[k]=ad[i]-*sigma*adb[i];
    }
  }

  method=2;

  DPSLDLT_Ordering(token,method);
  DPSLDLT_Preprocess(token,*neq,pointers,irowsgi,&ndim,&ops);

  if(*neq>200000){
    printf(" The out of core solver is used\n\n");
    DPSLDLT_OOCLimit(token,ooclimit);
    DPSLDLT_OOCPath(token,oocpath);
    DPSLDLT_FactorOOC(token,*neq,pointers,irowsgi,ausgi);
  }
  else{
    DPSLDLT_Factor(token,*neq,pointers,irowsgi,ausgi);
  }

  SFREE(pointers);

  return;
}

void sgi_solve(double *b,ITG token){

  DPSLDLT_Solve(token,b,b);

  return;
}

void sgi_cleanup(ITG token){

  DPSLDLT_Destroy(token);
  SFREE(irowsgi);
  SFREE(ausgi);

  return;
}

void sgi_main(double *ad, double *au, double *adb, double *aub, double *sigma,
         double *b, ITG *icol, ITG *irow, 
         ITG *neq, ITG *nzs, ITG token){

  if(*neq==0) return;

  sgi_factor(ad,au,adb,aub,sigma,icol,irow, 
             neq,nzs,token);

  sgi_solve(b,token);

  sgi_cleanup(token);

  return;
}

#endif

