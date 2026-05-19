!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
!     completing hel:
!
!     at the start of the subroutine: rhs of conservation of momentum
!                                     without pressure contribution
!     at the end of the subroutine: neighboring velocity terms subtracted
!
      subroutine complete_hel_blk(vel,hel,auv6,ipnei,neiel,nef,
     &        nactdohinv)
!
      implicit none
!
      integer i,j,k,indexf,ipnei(*),neiel(*),iel,nef,nactdohinv(*)
!
      real*8 hel(3,*),vel(nef,0:7),auv6(6,*)
!
      do i=1,nef
         indexf=ipnei(i)
         do j=1,6
            indexf=indexf+1
            iel=neiel(indexf)
            if(iel.eq.0) cycle
            do k=1,3
               hel(k,i)=hel(k,i)-auv6(j,i)*vel(iel,k)
c               write(*,*) 'complete_hel_blk ',nactdohinv(i),
c     &                nactdohinv(iel),k,
c     &                -auv6(j,i)
            enddo
         enddo
      enddo
!
c      do j=1,nef
c         write(*,*) 'complete_hel ',j,(hel(k,j),k=1,3)
c      enddo
!
      return
      end
