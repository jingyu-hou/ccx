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


void prediction(double *uam, ITG *nmethod, double *bet, double *gam, 
               double *dtime,
               ITG *ithermal, ITG *nk, double *veold, double *accold, double *v,
	       ITG *iinc, ITG *idiscon, double *vold, ITG *nactdof, ITG *mi){

    ITG j,k,mt=mi[1]+1;
    double dextrapol,scal1,scal2;
    FILE *pFile ;

    uam[0]=0.;
    uam[1]=0.;
    if(*nmethod==4){
      
	scal1=0.5*(1.-2.**bet)**dtime**dtime;
	scal2=(1.-*gam)**dtime;
	
	if(*ithermal<2){
	    for(k=0;k<*nk;++k){
		for(j=0;j<mt;j++){
		    dextrapol=*dtime*veold[mt*k+j]+scal1*accold[mt*k+j];
		    if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {uam[0]=fabs(dextrapol);}
		    v[mt*k+j]=vold[mt*k+j]+dextrapol;
		    veold[mt*k+j]=veold[mt*k+j]+scal2*accold[mt*k+j];
		    accold[mt*k+j]=0.;
		}
	    }
	}else if(*ithermal==2){
	    for(k=0;k<*nk;++k){
		for(j=0;j<mt;j++){
		    v[mt*k+j]=vold[mt*k+j];
		}
	    }
	    for(k=0;k<*nk;++k){
		dextrapol=*dtime*veold[mt*k];
		if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
		if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {uam[1]=fabs(dextrapol);}
		v[mt*k]+=dextrapol;
	    }
	}else{
	    for(k=0;k<*nk;++k){
		for(j=0;j<mt;++j){
		    dextrapol=*dtime*veold[mt*k+j]+scal1*accold[mt*k+j];
		    if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
		    if(j==0){
			if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {uam[1]=fabs(dextrapol);}
		    }else{
			if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {uam[0]=fabs(dextrapol);}
		    }
//		    v[mt*k+j]=vold[mt*k+j]+dextrapol;
//		    veold[mt*k+j]=veold[mt*k+j]+scal2*accold[mt*k+j];
		    v[mt*k+j]=vold[mt*k+j];
		    veold[mt*k+j]=veold[mt*k+j];
		    accold[mt*k+j]=0.;
		}
	    }
	}
    }
    
    /* for the static case: extrapolation of the previous increment
       (if any within the same step) */
    
    else{
	if(*iinc>1){
	    if(*ithermal<2){
		for(k=0;k<*nk;++k){
		    for(j=0;j<mt;++j){
			if(*idiscon==0){
			    dextrapol=*dtime*veold[mt*k+j];
			    if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {uam[0]=fabs(dextrapol);}	
			    v[mt*k+j]=vold[mt*k+j]+dextrapol;
			}else{
			    v[mt*k+j]=vold[mt*k+j];
			}
		    }
		}
	    }else if(*ithermal==2){
		for(k=0;k<*nk;++k){
		    for(j=0;j<mt;++j){
			v[mt*k+j]=vold[mt*k+j];
		    }
		}
		for(k=0;k<*nk;++k){
		    if(*idiscon==0){
			dextrapol=*dtime*veold[mt*k];
			if(fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
			if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k]>0)) {uam[1]=fabs(dextrapol);}	
			v[mt*k]+=dextrapol;
		    }
		}
	    }else{
		for(k=0;k<*nk;++k){
		    for(j=0;j<mt;++j){
			if(*idiscon==0){
			    dextrapol=*dtime*veold[mt*k+j];
			    if((j==0)&&fabs(dextrapol)>100.) dextrapol=100.*dextrapol/fabs(dextrapol);
			    if(j==0){
				if((fabs(dextrapol)>uam[1])&&(nactdof[mt*k+j]>0)) {uam[1]=fabs(dextrapol);}
			    }else{
				if((fabs(dextrapol)>uam[0])&&(nactdof[mt*k+j]>0)) {uam[0]=fabs(dextrapol);}
			    }	
			    v[mt*k+j]=vold[mt*k+j]+dextrapol;
			}else{
			    v[mt*k+j]=vold[mt*k+j];
			}
		    }
		}
	    }
	}
	else{
	    for(k=0;k<*nk;++k){
		for(j=0;j<mt;++j){
		    v[mt*k+j]=vold[mt*k+j];
		}
	    }
	}
    }
    *idiscon=0;

  return;
}
