/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdio.h>
#include <stdlib.h>
#include "WeICME.h"

int log_realloc=-1;

/*
 Diehl program
*/

void *u_calloc(size_t num,size_t size,const char *file,const int line, const char* ptr_name){

    /* allocating num elements of size bytes and initializing them to zero */

  void *a;
  char *env;

  if(num==0){
    a=NULL;
    return(a);
  }
      
  a=calloc(num,size);
  if(a==NULL){
    printf("*ERROR in u_calloc: error allocating memory\n");
    printf("variable=%s, file=%s, line=%d, num=%ld, size=%ld\n",ptr_name,file,line,num,size);
    if(num<0){
	printf("\n It looks like you may need the i8 (integer*8) version of WeICME\n");
    }
    exit(16);
  }
  else {
    if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
    }      
    if(log_realloc==1) {
	printf("ALLOCATION of variable %s, file %s, line=%d, num=%ld, size=%ld, address= %ld\n",ptr_name,file,line,num,size,(long int)a);
    }      
    return(a);
  }
}
