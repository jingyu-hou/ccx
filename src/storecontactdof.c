/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
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
#ifdef PARDISO
   #include "pardiso.h"
#endif

void storecontactdof(ITG *nope,ITG *nactdof, ITG *mt, ITG *konl, ITG **ikactcontp, 
          ITG *nactcont,ITG *nactcont_, double *bcont, double *fnl, 
          ITG *ikmpc, ITG *nmpc, ITG *ilmpc,ITG *ipompc, ITG *nodempc, 
          double *coefmpc){

  ITG j,j1,jdof,id,k,l,ist,index,node,ndir,*ikactcont=*ikactcontp;

  for(j=0;j<*nope;j++){
      for(j1=0;j1<3;j1++){
	  jdof=nactdof[*mt*(konl[j]-1)+j1+1];
	  if(jdof>0){
	      
	      jdof--;
	      FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
	      do{
		  if(id>0){
		      if(ikactcont[id-1]==jdof){
			  break;
		      }
		  }
		  (*nactcont)++;
		  if(*nactcont>*nactcont_){
		      *nactcont_=(ITG)(1.1**nactcont_);
		      RENEW(ikactcont,ITG,*nactcont_);
		  }
		  k=*nactcont-1;
		  l=k-1;
		  while(k>id){
		      ikactcont[k--]=ikactcont[l--];
		  }
		  ikactcont[id]=jdof;
		  break;
	      }while(1);
	      
	      bcont[jdof]-=fnl[3*j+j1];
	  }else{
	      jdof=8*(konl[j]-1)+j1+1;
	      FORTRAN(nident,(ikmpc,&jdof,nmpc,&id));
	      if(id>0){
		  if(ikmpc[id-1]==jdof){
		      id=ilmpc[id-1];
		      ist=ipompc[id-1];
		      index=nodempc[3*ist-1];
		      if(index==0) continue;
		      do{
			  node=nodempc[3*index-3];
			  ndir=nodempc[3*index-2];
			  jdof=nactdof[*mt*(node-1)+ndir];
			  if(jdof>0){
			      
			      jdof--;
			      FORTRAN(nident,(ikactcont,&jdof,nactcont,&id));
			      do{
				  if(id>0){
				      if(ikactcont[id-1]==jdof){
					  break;
				      }
				  }
				  (*nactcont)++;
				  if(*nactcont>*nactcont_){
				      *nactcont_=(ITG)(1.1**nactcont_);
				      RENEW(ikactcont,ITG,*nactcont_);
				  }
				  k=*nactcont-1;
				  l=k-1;
				  while(k>id){
				      ikactcont[k--]=ikactcont[l--];
				  }
				  ikactcont[id]=jdof;
				  break;
			      }while(1);
			      
/*			      bcont[jdof]+=coefmpc[index-1]*
			      fnl[3*j+j1]/coefmpc[ist-1];*/
			      bcont[jdof]-=coefmpc[index-1]*
				  fnl[3*j+j1]/coefmpc[ist-1];
			  }
			  index=nodempc[3*index-1];
			  if(index==0) break;
		      }while(1);
		  }
	      }
	  }
      }
  }

  *ikactcontp=ikactcont;
  
  return;
}

