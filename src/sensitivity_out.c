/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "WeICME.h"

void sensitivity_out(char *jobnamec,double *dgdxglob,ITG *neq,ITG *nobject,
       double *g0){
            
  char sensitivities[132]="",nominal[132]="";
  
  ITG i=0,iobject=0;
      
  FILE *f1;
           		 
  /* writing the sensitivities in the sen-file for optimizer */
        	
  strcpy(sensitivities,jobnamec);
  strcat(sensitivities,".sen");
  
  if((f1=fopen(sensitivities,"w"))==NULL){
      printf("*ERROR in sensitivity: cannot open sensitivity vector file for writing...");
      
      exit(0);
  }
  
  /* storing the sensitivity vectors */

  fprintf(f1,"---------------------------------- \n");
  fprintf(f1,"Objective \n");
  fprintf(f1,"---------------------------------- \n");
  
  for(i=0;i<neq[1];i++){
     for(iobject=0;iobject<*nobject;iobject++){
        fprintf(f1,"%12.5E",(double)dgdxglob[3+5*i+(5*neq[1]+2)*iobject]);
	fprintf(f1,";  ");    
     }
     fprintf(f1,"\n"); 
  }
  
  fclose(f1);

  /* writing the nominal values in the nom-file for optimizer */
        	
  strcpy(nominal,jobnamec);
  strcat(nominal,".nom");
  
  if((f1=fopen(nominal,"w"))==NULL){
      printf("*ERROR in sensitivity: cannot open sensitivity vector file for writing...");
      
      exit(0);
  }
  
  /* storing the sensitivity vectors */

  fprintf(f1,"---------------------------------- \n");
  fprintf(f1,"Objective \n");
  fprintf(f1,"---------------------------------- \n");
  
  for(iobject=0;iobject<*nobject;iobject++){
     fprintf(f1,"%12.5E",(double)g0[iobject]);
     fprintf(f1,";  ");    
  }
  fprintf(f1,"\n"); 
  
  fclose(f1);
  
  return;

}
