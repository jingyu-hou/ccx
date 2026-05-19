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

static ITG *nface1,*ielfa1,*ipnei1,*nef1,*num_cpus1;

static double *vel1,*flux1,*vfa1;

void hrk_udmain(ITG *nface,ITG *ielfa,double *vel,ITG *ipnei,ITG *nef,
		double *flux,double *vfa,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */
    
    ielfa1=ielfa;vel1=vel;ipnei1=ipnei;nef1=nef;flux1=flux;vfa1=vfa;
    num_cpus1=num_cpus;nface1=nface;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)hrk_ud1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of calcgammav1 */

void *hrk_ud1mt(ITG *i){

    ITG nfacea,nfaceb,nfacedelta;

    nfacedelta=(ITG)floor(*nface1/(double)(*num_cpus1));
    nfacea=*i*nfacedelta+1;
    nfaceb=(*i+1)*nfacedelta;
    if((*i==*num_cpus1-1)&&(nfaceb<*nface1)) nfaceb=*nface1;

    FORTRAN(hrk_ud,(ielfa1,vel1,ipnei1,nef1,flux1,vfa1,
			       &nfacea,&nfaceb));

    return NULL;
}
