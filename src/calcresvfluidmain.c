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

static ITG *ia1,*ja1,*nestart1,*nef1;

static double *a1,*b1,*au1,*x1,*res1=NULL,*xmax1=NULL;

void calcresvfluidmain(ITG *n,double *a,double *b,double *au,ITG *ia,
		       ITG *ja,double *x,double *res,ITG *nestart,ITG *num_cpus,
                       ITG *ncfd){

    ITG i;

    double xmax;
      
    /* variables for multithreading procedure */
    
    ITG *ithread=NULL;
    
    pthread_t tid[*num_cpus];

    /* residue of momentum equations in x */

    NNEW(res1,double,*num_cpus);

    a1=a;b1=b;au1=au;ia1=ia;ja1=ja;x1=x;nestart1=nestart;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcresvfluid1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    

    *res=0;
    for(i=0;i<*num_cpus;i++){
	*res+=res1[i];
    }
    SFREE(ithread);SFREE(res1);

    /* residue of momentum equations in y */

    NNEW(res1,double,*num_cpus);

    a1=a;b1=&b[*n];au1=au;ia1=ia;ja1=ja;x1=&x[*n];nestart1=nestart;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcresvfluid1mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    for(i=0;i<*num_cpus;i++){
	*res+=res1[i];
    }
    SFREE(ithread);SFREE(res1);

    /* residue of momentum equations in z */

    if(*ncfd>2){

	NNEW(res1,double,*num_cpus);

	a1=a;b1=&b[2**n];au1=au;ia1=ia;ja1=ja;x1=&x[2**n];nestart1=nestart;
	
	/* create threads and wait */
	
	NNEW(ithread,ITG,*num_cpus);
	for(i=0; i<*num_cpus; i++)  {
	    ithread[i]=i;
	    pthread_create(&tid[i], NULL, (void *)calcresvfluid1mt, (void *)&ithread[i]);
	}
	for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
	
	for(i=0;i<*num_cpus;i++){
	    *res+=res1[i];
	}
	SFREE(ithread);SFREE(res1);
    }

    /* max velocity */

    NNEW(xmax1,double,*num_cpus);

    x1=x;nef1=n;    
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)calcresvfluid2mt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
    
    xmax=0.;
    for(i=0;i<*num_cpus;i++){
	if(xmax1[i]>xmax){xmax=xmax1[i];}
    }
    SFREE(ithread);SFREE(xmax1);

    /* normalizing the momentum residue */

    *res=sqrt(*res/(*n))/xmax;
  
    return;

}

/* subroutine for multithreading of calcresvfluid1 */

void *calcresvfluid1mt(ITG *i){

    ITG n;

    /* number of equations for this thread */

    n=nestart1[*i+1]-nestart1[*i];

    FORTRAN(calcresvfluid1,(&n,a1,&b1[nestart1[*i]],&au1[nestart1[*i]],
			    ia1,&ja1[nestart1[*i]],&x1[nestart1[*i]],
			    &res1[*i]));

    return NULL;
}

/* subroutine for multithreading of calcresvfluid2 */

void *calcresvfluid2mt(ITG *i){

    ITG n;

    /* number of equations for this thread */

    n=nestart1[*i+1]-nestart1[*i];

    FORTRAN(calcresvfluid2,(&n,&x1[nestart1[*i]],&xmax1[*i],nef1));

    return NULL;
}
