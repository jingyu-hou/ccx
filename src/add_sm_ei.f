!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine add_sm_ei(au,ad,aub,adb,jq,irow,i,j,value,valuem,
     &  i0,i1)
!
!     stores the stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format and 
!     the mass coefficient (i,j) with value "valuem" in the lumped 
!     mass matrix
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1
!
      real*8 ad(*),au(*),adb(*),aub(*),value,valuem
!
      intent(in) jq,irow,i,j,value,valuem,
     &  i0,i1
!
      intent(inout) ad,au,adb,aub
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            ad(i)=ad(i)+value
            adb(i)=adb(i)+valuem
         else
            ad(i)=ad(i)+2.d0*value
            adb(i)=adb(i)+2.d0*valuem
         endif
         return
      elseif(i.gt.j) then
         ii=i
         jj=j
      else
         ii=j
         jj=i
      endif
c      write(*,*) ii,jj,value,valuem
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_ei: coefficient should be 0'
c         write(*,*) i,j,ii,jj,ipointer,irow(ipointer)
c         call exit(201)
      else
         au(ipointer)=au(ipointer)+value
         aub(ipointer)=aub(ipointer)+valuem
      endif
!
      return
      end













