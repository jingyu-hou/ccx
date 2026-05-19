!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine correctvel_simplec(adv,nef,volume,gradpel,vel,
     &  ipnei,auv)
!
!     correction of the velocity at the element centers due to the
!     pressure change (balance of mass)
!
!     the solution is stored in field bv.
!
      implicit none
!
      integer i,k,nef,ipnei(*),indexf
!
      real*8 adv(*),volume(*),gradpel(3,*),vel(nef,0:7),auv(*),a1
!
c$omp parallel default(none)
c$omp& shared(nef,vel,volume,gradpel,adv,ipnei,auv)
c$omp& private(i,k,indexf,a1)
c$omp do
      do i=1,nef
!
         a1=adv(i)
         do indexf=ipnei(i)+1,ipnei(i+1)
            a1=a1+auv(indexf)
         enddo
!
         do k=1,3
            vel(i,k)=vel(i,k)-volume(i)*gradpel(k,i)/a1
         enddo
      enddo
c$omp end do
c$omp end parallel
!  
      return
      end
