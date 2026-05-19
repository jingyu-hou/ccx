!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine orthotropic(orthol,anisox)
!
!     expands the 9 orthotropic elastic constants into a
!     3x3x3x3 matrix
!
      implicit none
!
      real*8 orthol(9),anisox(3,3,3,3)
!
      intent(in) orthol
!
      intent(out) anisox
!
      anisox(1,1,1,1)=orthol(1)
      anisox(1,1,1,2)=0.d0
      anisox(1,1,1,3)=0.d0
      anisox(1,1,2,1)=0.d0
      anisox(1,1,2,2)=orthol(2)
      anisox(1,1,2,3)=0.d0
      anisox(1,1,3,1)=0.d0
      anisox(1,1,3,2)=0.d0
      anisox(1,1,3,3)=orthol(4)
      anisox(1,2,1,1)=0.d0
      anisox(1,2,1,2)=orthol(7)
      anisox(1,2,1,3)=0.d0
      anisox(1,2,2,1)=orthol(7)
      anisox(1,2,2,2)=0.d0
      anisox(1,2,2,3)=0.d0
      anisox(1,2,3,1)=0.d0
      anisox(1,2,3,2)=0.d0
      anisox(1,2,3,3)=0.d0
      anisox(1,3,1,1)=0.d0
      anisox(1,3,1,2)=0.d0
      anisox(1,3,1,3)=orthol(8)
      anisox(1,3,2,1)=0.d0
      anisox(1,3,2,2)=0.d0
      anisox(1,3,2,3)=0.d0
      anisox(1,3,3,1)=orthol(8)
      anisox(1,3,3,2)=0.d0
      anisox(1,3,3,3)=0.d0
      anisox(2,1,1,1)=0.d0
      anisox(2,1,1,2)=orthol(7)
      anisox(2,1,1,3)=0.d0
      anisox(2,1,2,1)=orthol(7)
      anisox(2,1,2,2)=0.d0
      anisox(2,1,2,3)=0.d0
      anisox(2,1,3,1)=0.d0
      anisox(2,1,3,2)=0.d0
      anisox(2,1,3,3)=0.d0
      anisox(2,2,1,1)=orthol(2)
      anisox(2,2,1,2)=0.d0
      anisox(2,2,1,3)=0.d0
      anisox(2,2,2,1)=0.d0
      anisox(2,2,2,2)=orthol(3)
      anisox(2,2,2,3)=0.d0
      anisox(2,2,3,1)=0.d0
      anisox(2,2,3,2)=0.d0
      anisox(2,2,3,3)=orthol(5)
      anisox(2,3,1,1)=0.d0
      anisox(2,3,1,2)=0.d0
      anisox(2,3,1,3)=0.d0
      anisox(2,3,2,1)=0.d0
      anisox(2,3,2,2)=0.d0
      anisox(2,3,2,3)=orthol(9)
      anisox(2,3,3,1)=0.d0
      anisox(2,3,3,2)=orthol(9)
      anisox(2,3,3,3)=0.d0
      anisox(3,1,1,1)=0.d0
      anisox(3,1,1,2)=0.d0
      anisox(3,1,1,3)=orthol(8)
      anisox(3,1,2,1)=0.d0
      anisox(3,1,2,2)=0.d0
      anisox(3,1,2,3)=0.d0
      anisox(3,1,3,1)=orthol(8)
      anisox(3,1,3,2)=0.d0
      anisox(3,1,3,3)=0.d0
      anisox(3,2,1,1)=0.d0
      anisox(3,2,1,2)=0.d0
      anisox(3,2,1,3)=0.d0
      anisox(3,2,2,1)=0.d0
      anisox(3,2,2,2)=0.d0
      anisox(3,2,2,3)=orthol(9)
      anisox(3,2,3,1)=0.d0
      anisox(3,2,3,2)=orthol(9)
      anisox(3,2,3,3)=0.d0
      anisox(3,3,1,1)=orthol(4)
      anisox(3,3,1,2)=0.d0
      anisox(3,3,1,3)=0.d0
      anisox(3,3,2,1)=0.d0
      anisox(3,3,2,2)=orthol(5)
      anisox(3,3,2,3)=0.d0
      anisox(3,3,3,1)=0.d0
      anisox(3,3,3,2)=0.d0
      anisox(3,3,3,3)=orthol(6)
!
      return
      end

