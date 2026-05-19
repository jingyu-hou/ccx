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

static ITG *ipnei1,*nef1,*neifa1,*ielfa1,*neiel1,*ifabou1,*num_cpus1;

static double *area1,*vfa1,*xxna1,*flux1,*xxj1,*gradpfa1,*xlet1,*xle1,
    *vel1,*advfa1,*hfa1;

void calcfluxmain(double *area,double *vfa,double *xxna,ITG *ipnei,
		  ITG *nef,ITG *neifa,double *flux,double *xxj,
		  double *gradpfa,double *xlet,double *xle,double *vel,
		  double *advfa,ITG *ielfa,ITG *neiel,ITG *ifabou,
		  double *hfa,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */

    area1=area;vfa1=vfa;xxna1=xxna;ipnei1=ipnei;nef1=nef;neifa1=neifa;
    flux1=flux;xxj1=xxj;gradpfa1=gradpfa;xlet1=xlet;xle1=xle;vel1=vel;
    advfa1=advfa;ielfa1=ielfa;neiel1=neiel;ifabou1=ifabou;hfa1=hfa;
    num_cpus1=num_cpus;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcflux1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of materialdatacfd1 */

void *calcflux1mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(calcflux,(area1,vfa1,xxna1,ipnei1,nef1,neifa1,flux1,xxj1,
		       gradpfa1,xlet1,xle1,vel1,advfa1,ielfa1,neiel1,
		       ifabou1,hfa1,&nefa,&nefb));

    return NULL;
}
