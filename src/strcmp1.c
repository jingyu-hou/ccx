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

ITG strcmp1(const char *s1, const char *s2)
{
  ITG a,b;

  do {
    a=*s1++;
    b=*s2++;

/* the statement if((a=='\0')||(b=='\0')) has been treated separately
   in order to avoid the first field (s1) to be defined one longer
   than required; s1 is assumed to be a variable field, s2 is
   assumed to be a fixed string */

    if(b=='\0'){
      a='\0';
      b='\0';
      break;
    }
    if(a=='\0'){
      a='\0';
      b='\0';
      break;
    }
  }while(a==b);
  return(a-b);
}
	  
