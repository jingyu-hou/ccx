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

void writeheading(char *jobnamec,char *heading,ITG *nheading_){
    
    /* writes the headers in the frd-file */
    
    FILE *f1;

    char p1[6]="    1",fneig[132]="",c[2]="C",
         text[67]="                                                                  ";

    ITG i;
    
    strcpy(fneig,jobnamec);
    strcat(fneig,".frd");
    
    if((f1=fopen(fneig,"ab"))==NULL){
	printf("*ERROR in frd: cannot open frd file for writing...");
	exit(0);
    }

    /* first line */

    fprintf(f1,"%5s%1s\n",p1,c);

    /* header lines */

    for(i=0;i<*nheading_;i++){
	strcpy1(text,&heading[66*i],66);
	fprintf(f1,"%5sU%66s\n",p1,text);
    }

    fclose(f1);
    
    return;
}
