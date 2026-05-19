!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine add_bo_st(au,jq,irow,i,j,value)
!
!     stores the boundary stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format
!
      implicit none
!
      integer jq(*),irow(*),i,j,ipointer,id
!
      real*8 au(*),value
!
      intent(in) jq,irow,i,j,value
!
      intent(inout) au
!
      call nident(irow(jq(j)),i,jq(j+1)-jq(j),id)
!
      ipointer=jq(j)+id-1
!
      if(irow(ipointer).ne.i) then
c         write(*,*) i,j,ipointer,irow(ipointer)
         write(*,*) '*ERROR in add_bo_st: coefficient should be 0'
         call exit(201)
      else
         au(ipointer)=au(ipointer)+value
      endif
!
      return
      end













