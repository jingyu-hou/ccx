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

#include "readfrd.h"


/* liefert puffer aus string von position a bis b */
void stos(char *string, int a, int b, char *puffer)
{
  register int n, i;

  n=-1;
  for (i=a-1; i<b; i++)
    {
    n++;
    if ((i>=MAX_LINE_LENGTH)||(n>=MAX_LINE_LENGTH)) break;
    puffer[n] = string[i];
  }
  puffer[n+1] = '\0';
}

/* schreibt string in puffer von position a bis b */
void stos_inv(char *string, int a, int b, char *puffer)
{
  register int n, i;

  n=-1;
  for (i=a-1; i<b; i++)
    {
    n++;
    if ((i>132)||(n>132)) break;
    puffer[i] = string[n];
  }
}

