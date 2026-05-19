!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine create_iau6(nef,ipnei,neiel,jq,irow,nzs,iau6,lakonf)
!
!     sets up a field of pointers iau6(j,i) for neighbor j of
!     element (cell) i (for CFD-applications) into field auv,aup,aut
!
      implicit none
!
      character*8 lakonf(*)
!
      integer nef,ipnei(*),neiel(*),jq(*),irow(*),nzs,iau6(6,*),
     &  numfaces,id,i,j,iel,indexf
!
!
!
      do i=1,nef
         indexf=ipnei(i)
!
         do j=1,ipnei(i+1)-ipnei(i)
            indexf=indexf+1
            iel=neiel(indexf)
            if(iel.eq.0) cycle
            if(i.gt.iel) then
               call nident(irow(jq(iel)),i,jq(iel+1)-jq(iel),id)
               iau6(j,i)=jq(iel)+id-1
            elseif(i.lt.iel) then
               call nident(irow(jq(i)),iel,jq(i+1)-jq(i),id)
               iau6(j,i)=nzs+jq(i)+id-1
            endif
         enddo
      enddo
!     
      return
      end
