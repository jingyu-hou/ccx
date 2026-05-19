/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdio.h>
#include <stdlib.h>
#include "WeICME.h"
extern int log_realloc;

/*
 Diehl program
*/

void *u_realloc(void* ptr,size_t size,const char *file,const int line, const char* ptr_name){

    /* reallocating a field with pointer ptr to size bytes */

  void *a;
  char *env;

  a=realloc(ptr,size);

  if(a==NULL && ptr!=NULL && size!=0){
    printf("*ERROR in u_realloc: error allocating memory\n");
    printf("variable=%s, file=%s, line=%d, size(bytes)=%ld, oldaddress=%ld\n",ptr_name,file,line,size,(long int)ptr);
    exit(16);
  }
  else {
    if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
    }      
    if(log_realloc==1) {
      printf("REALLOCATION of variable %s, file %s, line=%d: size(bytes)=%ld, oldaddress= %ld,address= %ld\n",ptr_name,file,line,size,(long int)ptr,(long int)a);
    }      
    return(a);
  }
}
