!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!  
      subroutine add_sm_st_corio(au,ad,jq,irow,i,j,value,i0,i1)
!
!     stores the stiffness coefficient (i,j) with value "value"
!     in the stiffness matrix stored in spare matrix format
!
!     modification for Coriolis: the Coriolis matrix is antisymmetric,
!     i.e. the transpose is the negative matrix: A^T=-A
!
      implicit none
!
      integer jq(*),irow(*),i,j,ii,jj,ipointer,id,i0,i1
      real*8 ad(*),au(*),value,valuenew
!
      if(i.eq.j) then
c         if(i0.eq.i1) then
c            ad(i)=ad(i)+value
c         else
c            ad(i)=ad(i)+2.d0*value
c         endif
         return
      elseif(i.gt.j) then
         ii=i
         jj=j
         valuenew=value
      else
         ii=j
         jj=i
         valuenew=-value
      endif
!
      call nident(irow(jq(jj)),ii,jq(jj+1)-jq(jj),id)
!
      ipointer=jq(jj)+id-1
!
      if(irow(ipointer).ne.ii) then
         write(*,*) '*ERROR in add_sm_st: coefficient should be 0'
         call exit(201)
      else
         au(ipointer)=au(ipointer)+valuenew
      endif
!
      return
      end













