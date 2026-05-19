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

static ITG *iamorig1,*nz_num1,*num_cpus1;

static double *au1,*am1;

void reorderlhsmain(double *au,double *am,ITG *iamorig,ITG *nz_num,
		    ITG *num_cpus){

    ITG i;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;;
    
    pthread_t tid[*num_cpus];

    /* calculation of the density at the cell centers */
    
    au1=au;am1=am;iamorig1=iamorig;nz_num1=nz_num;num_cpus1=num_cpus;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)reorderlhs1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);
  
  return;

}

/* subroutine for multithreading of calcgammav1 */

void *reorderlhs1mt(ITG *i){

    ITG nz_numa,nz_numb,nz_numdelta;

    nz_numdelta=(ITG)floor(*nz_num1/(double)(*num_cpus1));
    nz_numa=*i*nz_numdelta+1;
    nz_numb=(*i+1)*nz_numdelta;
    if((*i==*num_cpus1-1)&&(nz_numb<*nz_num1)) nz_numb=*nz_num1;

    FORTRAN(reorderlhs,(au1,am1,iamorig1,&nz_numa,&nz_numb));

    return NULL;
}
