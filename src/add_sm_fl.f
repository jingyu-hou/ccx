!
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine add_sm_fl(aub,adb,jq,irow,i,j,value,
     &  i0,i1)
!
!     stores the coefficient (i,j) with value "value" in the
!     fluid matrix
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1
      real*8 adb(*),aub(*),value
!
      if(i.eq.j) then
         if(i0.eq.i1) then
            adb(i)=adb(i)+value
         else
            adb(i)=adb(i)+2.d0*value
         endif
         return
      elseif(i.gt.j) then
         ii=i
         jj=j
      else
         ii=j
         jj=i
      endif
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_ei: coefficient should be 0'
c         write(*,*) i,j,ii,jj,ipointer,irow(ipointer)
      else
         aub(ipointer)=aub(ipointer)+value
      endif
!
      return
      end













