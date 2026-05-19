/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
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

static char *objectset1;

static ITG *nobject1,*nk1,*nodedesi1,*ndesi1,*nx1,*ny1,*nz1,*neighbor1=NULL,
    num_cpus;

/* y1 had to be replaced by yy1, else the following compiler error
   popped up: 

   filtermain.c:42: error: ‘y1’ redeclared as different kind of symbol */


static double *dgdxglob1,*xo1,*yo1,*zo1,*x1,*yy1,*z1,*r1=NULL,*xdesi1,
              *distmin1;

void filtermain(double *co, double *dgdxglob, ITG *nobject, ITG *nk,
                ITG *nodedesi, ITG *ndesi, char *objectset,double *xdesi,
		double *distmin){

    /* filtering the sensitivities */

    ITG *nx=NULL,*ny=NULL,*nz=NULL,i,*ithread=NULL;
    
    double *xo=NULL,*yo=NULL,*zo=NULL,*x=NULL,*y=NULL,*z=NULL;
    
    /* if no radius is defined no filtering is performed
       the radius applies to all objective functions */
    
    if(*nobject==0){return;}
    if(strcmp1(&objectset[81],"     ")==0){
	for(i=1;i<2**nk**nobject;i=i+2){
	    dgdxglob[i]=dgdxglob[i-1];
	}
	return;
    }
    
    /* prepare for near3d_se */
    
    NNEW(xo,double,*ndesi);
    NNEW(yo,double,*ndesi);
    NNEW(zo,double,*ndesi);
    NNEW(x,double,*ndesi);
    NNEW(y,double,*ndesi);
    NNEW(z,double,*ndesi);
    NNEW(nx,ITG,*ndesi);
    NNEW(ny,ITG,*ndesi);
    NNEW(nz,ITG,*ndesi);
    
    FORTRAN(prefilter,(co,nodedesi,ndesi,xo,yo,zo,x,y,z,nx,ny,nz));
    
    /* variables for multithreading procedure */
    
    ITG sys_cpus;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;
    
    /* explicit user declaration prevails */
    
    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
	sys_cpus=atoi(envsys);
	if(sys_cpus<0) sys_cpus=0;
    }
    
    /* automatic detection of available number of processors */
    
    if(sys_cpus==0){
	sys_cpus = getSystemCPUs();
	if(sys_cpus<1) sys_cpus=1;
    }
    
    /* local declaration prevails, if strictly positive */
    
    envloc = getenv("CCX_NPROC_SENS");
    if(envloc){
	num_cpus=atoi(envloc);
	if(num_cpus<0){
	    num_cpus=0;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
	
    }
    
    /* else global declaration, if any, applies */
    
    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
	if (env)
	    num_cpus = atoi(env);
	if (num_cpus < 1) {
	    num_cpus=1;
	}else if(num_cpus>sys_cpus){
	    num_cpus=sys_cpus;
	}
    }
    
    /* check that the number of cpus does not supercede the number
       of design variables */
    
    if(*ndesi<num_cpus) num_cpus=*ndesi;
    
    pthread_t tid[num_cpus];
    
    NNEW(neighbor1,ITG,num_cpus*(*ndesi+6));
    NNEW(r1,double,num_cpus*(*ndesi+6));
    
    dgdxglob1=dgdxglob;nobject1=nobject;nk1=nk;nodedesi1=nodedesi;
    ndesi1=ndesi;objectset1=objectset;xo1=xo;yo1=yo;zo1=zo;
    x1=x;yy1=y;z1=z;nx1=nx;ny1=ny;nz1=nz;xdesi1=xdesi;
    distmin1=distmin;
    
    /* filtering */
    
    printf(" Using up to %" ITGFORMAT " cpu(s) for filtering the sensitivities.\n\n", num_cpus);
    
    /* create threads and wait */
  
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)filtermt, (void *)&ithread[i]);
    }
    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(neighbor1);SFREE(r1);SFREE(xo);SFREE(yo);SFREE(zo);
    SFREE(x);SFREE(y);SFREE(z);SFREE(nx);SFREE(ny);SFREE(nz);
    SFREE(ithread);
    
    return;
    
} 

/* subroutine for multithreading of filter */

void *filtermt(ITG *i){

    ITG indexr,ndesia,ndesib,ndesidelta;

    indexr=*i*(*ndesi1+6);
    
    ndesidelta=(ITG)ceil(*ndesi1/(double)num_cpus);
    ndesia=*i*ndesidelta+1;
    ndesib=(*i+1)*ndesidelta;
    if(ndesib>*ndesi1) ndesib=*ndesi1;

//    printf("i=%" ITGFORMAT ",ntria=%" ITGFORMAT ",ntrib=%" ITGFORMAT "\n",i,ntria,ntrib);
//    printf("indexad=%" ITGFORMAT ",indexau=%" ITGFORMAT ",indexdi=%" ITGFORMAT "\n",indexad,indexau,indexdi);

    FORTRAN(filter,(dgdxglob1,nobject1,nk1,nodedesi1,ndesi1,objectset1,
                    xo1,yo1,zo1,x1,yy1,z1,nx1,ny1,nz1,&neighbor1[indexr],
                    &r1[indexr],&ndesia,&ndesib,xdesi1,distmin1));

    return NULL;
}
    
