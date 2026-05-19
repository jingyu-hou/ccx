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

/*---------------------------------------------------------------------*/
/* Strings vergleichen (bis zu welchem character sind sie gleich?)     */
/*---------------------------------------------------------------------*/

int compare (char *str1, char *str2, int length)
{
    int     i;

    i = 0;
    while ((str1[i]==str2[i]) && (i<length))
        i++;

    return i;               /* Return how far we got before a difference occurred, or the variable      */
                            /* length, whichever is the smaller                                         */
}

