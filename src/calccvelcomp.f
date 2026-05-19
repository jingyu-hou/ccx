!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calccvelcomp(nef,vel,shcon,nshcon,ielmat,ntmat_,
     &  mi,cvel,physcon)
!
!     calculation of the secant heat capacity at constant volume
!     in the element centers (compressible media)
!
      implicit none
!
      integer nef,i,nshcon(*),imat,ntmat_,mi(*),ielmat(mi(3),*)
!
      real*8 t1l,vel(nef,0:7),cp,shcon(0:3,ntmat_,*),cvel(*),
     &  physcon(*)
!     
c$omp parallel default(none)
c$omp& shared(nef,vel,ielmat,ntmat_,shcon,nshcon,physcon,cvel)
c$omp& private(i,t1l,imat,cp)
c$omp do
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmat(1,i)
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvel(i)=cp-shcon(3,1,imat)
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
