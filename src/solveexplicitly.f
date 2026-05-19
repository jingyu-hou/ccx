!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine solveexplicitly(nef,vel,bv,auv,ipnei,neiel,nflnei)
!
!     explicitly solving the momentum equations
!     the solution is stored in the rhs
!
      implicit none
!
      integer i,j,indexf,ipnei(*),iel,neiel(*),nflnei,nef
!
      real*8 vel(nef,0:7),bv(nef,3),auv(*)
!
c$omp parallel default(none)
c$omp& shared(nef,vel,bv,neiel,nflnei,auv,ipnei)
c$omp& private(i,j,indexf,iel)
c$omp do
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
            iel=neiel(indexf)
            do j=1,3
               bv(i,j)=bv(i,j)-auv(indexf)*vel(iel,j)
            enddo
         enddo
         do j=1,3
            bv(i,j)=bv(i,j)/auv(nflnei+i)
         enddo
      enddo
c$omp end do
c$omp end parallel
! 
      return
      end
