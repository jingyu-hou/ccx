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

void insertcbs(ITG *ipointer, ITG **irowp, ITG **nextp, ITG *i1,
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

       notice that in C the positions start at 0 and not at 1 as in 
       FORTRAN; the present routine is written in FORTRAN convention 

       routine specifically for the CBS-Method in CFD */
  
  ITG idof1,idof2,istart,*irow=NULL,*next=NULL;

  irow=*irowp;
  next=*nextp;

  if(*i1<*i2){
    idof1=*i1;
    idof2=*i2;
  }
  else{
    idof1=*i2;
    idof2=*i1;
  }

  if(ipointer[idof2-1]==0){
    ++*ifree;
    if(*ifree>*nzs_){
      *nzs_=(ITG)(1.1**nzs_);
      RENEW(irow,ITG,*nzs_);
      RENEW(next,ITG,*nzs_);
    }
    ipointer[idof2-1]=*ifree;
    irow[*ifree-1]=idof1;
    next[*ifree-1]=0;
  }
  else{
    istart=ipointer[idof2-1];
    while(1){
      if(irow[istart-1]==idof1) break;
      if(next[istart-1]==0){
	++*ifree;
	if(*ifree>*nzs_){
	  *nzs_=(ITG)(1.1**nzs_);
	  RENEW(irow,ITG,*nzs_);
	  RENEW(next,ITG,*nzs_);
	}
	next[istart-1]=*ifree;
	irow[*ifree-1]=idof1;
	next[*ifree-1]=0;
	break;
      }
      else{
	istart=next[istart-1];
      }
    }
  }

  *irowp=irow;
  *nextp=next;
  
  return;

}
	  
/*

Here starts the original FORTRAN code, which was transferred to the
C-code above in order to allow automatic reallocation

      subroutine insert(ipointer,irow,next,i1,i2,ifree,nzs_)
!
!     inserts a new nonzero matrix position into the data structure
!
      implicit none
!
      integer ipointer(*),irow(*),next(*),i1,i2,ifree,nzs_,idof1,
     &  idof2,istart
!
      if(i1.lt.i2) then
        idof1=i1
        idof2=i2
      else
        idof1=i2
        idof2=i1
      endif
!
      if(ipointer(idof2).eq.0) then
        ifree=ifree+1
        if(ifree.gt.nzs_) then
           write(*,*) '*ERROR in insert: increase nzs_'
           stop
        endif
        ipointer(idof2)=ifree
        irow(ifree)=idof1
        next(ifree)=0
      else
        istart=ipointer(idof2)
        do
          if(irow(istart).eq.idof1) exit
          if(next(istart).eq.0) then
            ifree=ifree+1
            if(ifree.gt.nzs_) then
               write(*,*) '*ERROR in insert: increase nzs_'
               stop
            endif
            next(istart)=ifree
            irow(ifree)=idof1
            next(ifree)=0
            exit
          else
            istart=next(istart)
          endif
        enddo
      endif
!
      return
      end

      */
