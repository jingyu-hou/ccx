!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine add_sm_fl_as(au,ad,jq,irow,i,j,value,nzs)
!
!     stores the stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format
!     asymmetric version for fluid dynamics
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,nzs,ioffset
!
      real*8 ad(*),au(*),value
!
      intent(in) jq,irow,i,j,value,nzs
!
      intent(inout) au,ad
!
      if(i.eq.j) then
         ad(i)=ad(i)+value
         return
      elseif(i.gt.j) then
         ioffset=0
         ii=i
         jj=j
      else
         ioffset=nzs
         ii=j
         jj=i
      endif
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_fl_as: coefficient should be 0'
         call exit(201)
      else
         ipointer=ipointer+ioffset
c         write(*,*) 'add_sm_fl_as ',i,j,ipointer
         au(ipointer)=au(ipointer)+value
      endif
!
      return
      end













