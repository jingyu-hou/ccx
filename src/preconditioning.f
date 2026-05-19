!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine preconditioning(ad,au,b,neq,irow,jq,adaux)
!
      implicit none
!
      integer neq,irow(*),jq(*),i,ic,ir
!
      real*8 ad(*),au(*),b(*),adaux(*),adc
!
      intent(in) neq,irow,jq
!
      intent(inout) ad,au,b,adaux
!
!     inverse of the square root of the diagonal
!     the sign takes care that the diagonal term becomes 1 
!     (and not -1)
!
      do i=1,neq
         adaux(i)=dsign(1.d0/dsqrt(dabs(ad(i))),ad(i))
      enddo
c      do i=1,neq
c         adaux(i)=1.d0/dsqrt(dabs(ad(i)))
c      enddo
!
!     scaling the matrix and the right hand side
!
      do ic=1,neq
         adc=dabs(adaux(ic))
!
!        scaling the diagonal
!
         ad(ic)=ad(ic)*adc*adaux(ic)
!
!        scaling the off-diagonal terms
!
         do i=jq(ic),jq(ic+1)-1
            ir=irow(i)
            au(i)=au(i)*adc*adaux(ir)
         enddo
!
!        scaling the right hand side
!
         b(ic)=b(ic)*adaux(ic)
      enddo
!
!     taking the absolute value
!
      do i=1,neq
         adaux(i)=dabs(adaux(i))
      enddo
!
      return
      end
