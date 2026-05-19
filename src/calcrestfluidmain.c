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

static ITG *ia1,*ja1,*nestart1;

static double *a1,*b1,*au1,*x1,*res1=NULL,*xmin1=NULL,*xmax1=NULL;

void calcrestfluidmain(ITG *n,double *a,double *b,double *au,ITG *ia,
		ITG *ja,double *x,double *res,ITG *nestart,ITG *num_cpus){

    ITG i;

    double xmin,xmax;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;
    
    pthread_t tid[*num_cpus];

    /* residue of momentum equations in x */

    NNEW(res1,double,*num_cpus);
    NNEW(xmin1,double,*num_cpus);
    NNEW(xmax1,double,*num_cpus);

    a1=a;b1=b;au1=au;ia1=ia;ja1=ja;x1=x;nestart1=nestart;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcrestfluidmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    SFREE(ithread);

    /* total residual */

    *res=0;
    for(i=0;i<*num_cpus;i++){
	*res+=res1[i];
    }
    SFREE(res1);
    
    /* maximum and minimum temperature */

    xmax=0.;
    xmin=0.;  /* is this correct? not 1.e30 ? */
    for(i=0;i<*num_cpus;i++){
	if(xmax1[i]>xmax){xmax=xmax1[i];}
	if(xmin1[i]<xmin){xmin=xmin1[i];}
    }
    SFREE(xmax1);SFREE(xmin1);

    /* normalizing the momentum residue */

    *res=sqrt(*res/(*n))/(xmax-xmin);
  
    return;

}

/* subroutine for multithreading of calcresvfluid1 */

void *calcrestfluidmt(ITG *i){

    ITG n;

    /* number of equations for this thread */

    n=nestart1[*i+1]-nestart1[*i];

    FORTRAN(calcrestfluid,(&n,a1,&b1[nestart1[*i]],&au1[nestart1[*i]],
			   ia1,&ja1[nestart1[*i]],&x1[nestart1[*i]],
			   &res1[*i],&xmin1[*i],&xmax1[*i]));

    return NULL;
}
