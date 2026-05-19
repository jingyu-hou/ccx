!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine mafillpbc(nef,au,ad,jq,irow,
     &  b,iatleastonepressurebc,nzs)
!
!     filling the lhs and rhs to calculate p
!
      implicit none
!
      integer i,nef,irow(*),jq(*),iatleastonepressurebc,nzs
!
      real*8 ad(*),au(*),b(*)
!     
!     at least one pressure bc is needed. If none is applied,
!     the last dof is set to 0
!     
!     a pressure bc is only recognized if not all velocity degrees of
!     freedom are prescribed on the same face
!     
c      write(*,*) 'mafillpbc', iatleastonepressurebc
      if(iatleastonepressurebc.eq.0) then
         ad(nef)=1.d0
         b(nef)=0.d0
         do i=2,nef
            if(jq(i)-1>0) then
               if(irow(jq(i)-1).eq.nef) then
                  au(jq(i)-1)=0.d0
               endif
            endif
         enddo
      endif
!     
c      do i=1,nzs
c         write(*,*) 'mafillp irow,au',i,au(i)
c      enddo
c      do i=1,nef
c         write(*,*) 'mafillp ad b',i,ad(i),b(i)
c      enddo
!     
      return
      end
