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

void *u_free(void* ptr,const char *file,const int line, const char* ptr_name){

    /* freeing a field with pointer ptr  */

  char *env;

  free(ptr);

  if(log_realloc==-1) {
      log_realloc=0;
      env=getenv("CCX_LOG_ALLOC");
      if(env) {log_realloc=atoi(env);}
  }      
  if(log_realloc==1) {
      printf("FREEING of variable %s, file %s, line=%d: oldaddress= %ld\n",ptr_name,file,line,(long int)ptr);
  }      
  return;
}
