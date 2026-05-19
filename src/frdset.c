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

void frdset(char *filabl,char *set,ITG *iset,ITG *istartset,ITG *iendset,
	    ITG *ialset,ITG *inum,ITG *noutloc,ITG *nout,ITG *nset,
	    ITG *noutmin,ITG *noutplus,ITG *iselect,ITG *ngraph){

  ITG j,k;

  char noset[81];
     
  /* check for a set, if any */

  strcpy1(noset,&filabl[6],81);
  for((*iset)=0;(*iset)<(*nset);(*iset)++){
    if(strcmp2(&set[81**iset],noset,81)==0) break;
  }
  (*iset)++;
  if(*iset>*nset)*iset=0;
  //    printf("iset,noutplus %" ITGFORMAT " %" ITGFORMAT "\n",*iset,*noutplus);

  /* determining the number of nodes in the set */

  if(*iset==0){

    /* no set defined */

    //    printf("iselect,noutplus %" ITGFORMAT " %" ITGFORMAT "\n",*iselect,*noutplus);

    if(*iselect==1){
      *noutloc=*noutplus;
    }else if(*iselect==-1){
      *noutloc=*noutmin;
    }else{
      *noutloc=*nout;
    }

  }else{

    /* a set was defined */

    *noutloc=0;
    for(j=istartset[*iset-1]-1;j<iendset[*iset-1];j++){
      if(ialset[j]>0){
	if(*iselect==-1){
	  if(inum[ialset[j]-1]<0) (*noutloc)++;
	}else if(*iselect==1){
	  if(inum[ialset[j]-1]>0) (*noutloc)++;
	}else{
	  if(inum[ialset[j]-1]!=0) (*noutloc)++;
	}
      }else{
	k=ialset[j-2];
	do{
	  k=k-ialset[j];
	  if(k>=ialset[j-1]) break;
	  if(*iselect==-1){
	    if(inum[k-1]<0) (*noutloc)++;
	  }else if(*iselect==1){
	    if(inum[k-1]>0) (*noutloc)++;
	  }else{
	    if(inum[k-1]!=0) (*noutloc)++;
	  }
	}while(1);
      }
    }
    if(*ngraph>1) (*noutloc)*=(*ngraph);
  }
    

}


     /*!
!
!
      subroutine frdset(filabl,set,iset,istartset,iendset,
     &  ialset,inum,noutloc,nout,nset,noutmin,noutplus,iselect,
     &  ngraph)
!
!     stores the results in frd format
!
      implicit none
!
      character*81 set(*),noset
      character*87 filabl
!
      integer iset,istartset(*),iendset(*),ialset(*),inum(*),
     &  noutloc,j,k,nout,nset,noutmin,noutplus,iselect,ngraph
!     
!     check for a set, if any
!     
      noset=filabl(7:87)
      do iset=1,nset
         if(set(iset).eq.noset) exit
      enddo
      if(iset.gt.nset) iset=0
!     
!     determining the number of nodes in the set
!     
      if(iset.eq.0) then
         if(iselect.eq.1) then
            noutloc=noutplus
         elseif(iselect.eq.-1) then
            noutloc=noutmin
         else
            noutloc=nout
         endif
      else
         noutloc=0
         do j=istartset(iset),iendset(iset)
            if(ialset(j).gt.0) then
               if(iselect.eq.-1) then
                  if(inum(ialset(j)).lt.0) noutloc=noutloc+1
               elseif(iselect.eq.1) then
                  if(inum(ialset(j)).gt.0) noutloc=noutloc+1
               else
                  if(inum(ialset(j)).ne.0) noutloc=noutloc+1
               endif
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  if(iselect.eq.-1) then
                     if(inum(k).lt.0) noutloc=noutloc+1
                  elseif(iselect.eq.1) then
                     if(inum(k).gt.0) noutloc=noutloc+1
                  else
                     if(inum(k).ne.0) noutloc=noutloc+1
                  endif
               enddo
            endif
         enddo
         if(ngraph.gt.1) noutloc=noutloc*ngraph
      endif
!     
      return
      end*/


