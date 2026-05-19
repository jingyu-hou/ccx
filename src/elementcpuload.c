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

void elementcpuload(ITG *neapar,ITG *nebpar,ITG *ne,ITG *ipkon,ITG *num_cpus){

 /*  divides the elements into ranges with an equal number of
     active elements (element numbering may have gaps) for
     parallel processing on different cpus */

    ITG i,nepar,*ipar=NULL,idelta,isum;
    
    NNEW(ipar,ITG,*ne);

    nepar=0;
    for(i=0;i<*ne;i++){
	if(ipkon[i]>-1){

	    /* active element */
	    
	    ipar[nepar]=i;
	    nepar++;
	}
    }
    if(nepar<*num_cpus) *num_cpus=nepar;

    /* dividing the element number range into num_cpus equal numbers of 
       active elements */

    idelta=nepar/(*num_cpus);
    isum=0;
    for(i=0;i<*num_cpus;i++){
	neapar[i]=ipar[isum];
	if(i!=*num_cpus-1){
	    isum+=idelta;
	}else{
	    isum=nepar;
	}
	nebpar[i]=ipar[isum-1];
    }
    
    SFREE(ipar);
    
    return;
}
