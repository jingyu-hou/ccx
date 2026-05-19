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

void insertrad(ITG *ipointer, ITG **irowp, ITG **nextp, ITG *i1,
	    ITG *i2, ITG *ifree, ITG *nzs_){

  /*   inserts a new nonzero matrix position into the data structure 
       in FORTRAN notation: 
       - ipointer(i) points to a position in field irow containing
         the row number of a nonzero position in column i; 
         next(ipointer(i)) points a position in field irow containing
         the row number of another nonzero position in column i, and
         so on until no nonzero positions in column i are left; for 
         the position j in field irow containing the momentarily last
         nonzero number in column i we have next(j)=0 

         special version of insert.c for the call in mastructrad.c

       notice that in C the positions start at 0 and not at 1 as in 
       FORTRAN; the present routine is written in FORTRAN convention */

  ITG *irow=NULL,*next=NULL;

  irow=*irowp;
  next=*nextp;

  ++*ifree;
  if(*ifree>*nzs_){
      *nzs_=(ITG)(1.1**nzs_);
      RENEW(irow,ITG,*nzs_);
      RENEW(next,ITG,*nzs_);
  }
  
  irow[*ifree-1]=*i2;
  next[*ifree-1]=ipointer[*i1-1];
  ipointer[*i1-1]=*ifree;

  *irowp=irow;
  *nextp=next;
  
  return;

}
