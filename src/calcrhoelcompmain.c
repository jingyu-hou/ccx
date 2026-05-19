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

static ITG *nef1,*ielmatf1,*ntmat1_,*mi1,*num_cpus1;

static double *vel1,*shcon1;

void calcrhoelcompmain(ITG *nef,double *vel,double *shcon,ITG *ielmatf,
		       ITG *ntmat_,ITG *mi,ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */

    nef1=nef;vel1=vel;shcon1=shcon;ielmatf1=ielmatf;ntmat1_=ntmat_;mi1=mi;
    num_cpus1=num_cpus;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcrhoelcomp1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of materialdatacfd1 */

void *calcrhoelcomp1mt(ITG *i){

    ITG nefa,nefb,nefdelta;

    nefdelta=(ITG)floor(*nef1/(double)(*num_cpus1));
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
    if((*i==*num_cpus1-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(calcrhoelcomp,(nef1,vel1,shcon1,ielmatf1,ntmat1_,mi1,
			    &nefa,&nefb));

    return NULL;
}
