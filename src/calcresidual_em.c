/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdio.h>
#include <math.h> 
#include <stdlib.h>
#include "WeICME.h"
#ifdef SPOOLES
   #include "spooles.h"
#endif
#ifdef SGI
   #include "sgi.h"
#endif
#ifdef TAUCS
   #include "tau.h"
#endif


void calcresidual_em(ITG *nmethod, ITG *neq, double *b, double *fext, double *f,
        ITG *iexpl, ITG *nactdof, double *aux1, double *aux2, double *vold,
        double *vini, double *dtime, double *accold, ITG *nk, double *adb,
        double *aub, ITG *jq, ITG *irow, ITG *nzl, double *alpha,
        double *fextini, double *fini, ITG *islavnode, ITG *nslavnode,
        ITG *mortar, ITG *ntie,double *f_cm,
	double* f_cs, ITG *mi,ITG *nzs,ITG *nasym,ITG *ithermal){

    ITG j,k,mt=mi[1]+1,jstart;
      
    /* residual for a static analysis */
      
    if(*nmethod!=4){
	for(k=0;k<neq[1];++k){
	    b[k]=fext[k]-f[k];
	}
    }
      
    /* residual for implicit dynamics */
      
    else{

	if(*ithermal<2){
	    jstart=1;
	}else{
	    jstart=0;
	}

        /* calculating a pseudo-velocity */

/*	for(k=0;k<*nk;++k){
	    for(j=jstart;j<mt;++j){
		if(nactdof[mt*k+j]>0){aux2[nactdof[mt*k+j]-1]=(vold[mt*k+j]-vini[mt*k+j])/(*dtime);}
	    }
	    }*/

	for(k=0;k<*nk;++k){
	    for(j=jstart;j<1;++j){
		if(nactdof[mt*k+j]>0){aux2[nactdof[mt*k+j]-1]=(vold[mt*k+j]-vini[mt*k+j])/(*dtime);}
	    }
	}

	for(k=0;k<*nk;++k){
	    for(j=1;j<mt;++j){
		if(nactdof[mt*k+j]>0){aux2[nactdof[mt*k+j]-1]=(-vini[mt*k+j])/(*dtime);}
	    }
	}

        /* calculating "capacity"-matrix times pseudo-velocity */

	if(*nasym==0){
	    FORTRAN(op,(&neq[1],aux2,b,adb,aub,jq,irow)); 
	}else{
	    FORTRAN(opas,(&neq[1],aux2,b,adb,aub,jq,irow,nzs)); 
	}

	for(k=0;k<neq[1];++k){
	    b[k]=fext[k]-f[k]-b[k];
	} 
    }

    return;
}
