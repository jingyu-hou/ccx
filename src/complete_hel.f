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
      subroutine complete_hel(nef,bv,hel,adv,auv,ipnei,neiel,nzs)
!
      implicit none
!
      integer neiel(*),nef,nzs,j,k,l,jdof1,ipnei(*),indexf,i,iel
!
      real*8 hel(3,*),bv(nef,3),auv(*),adv(*)
!
!     off-diagonal terms
!
c$omp parallel default(none)
c$omp& shared(nef,ipnei,neiel,hel,auv,bv)
c$omp& private(i,indexf,iel,k)
c$omp do
      do i=1,nef
         do indexf=ipnei(i)+1,ipnei(i+1)
!
!           neighboring element
!
            iel=neiel(indexf)
            if(iel.ne.0) then
               do k=1,3
                  hel(k,i)=hel(k,i)-auv(indexf)*bv(iel,k)
               enddo
            endif
         enddo
      enddo
c$omp end do
c$omp end parallel
!
      return
      end
