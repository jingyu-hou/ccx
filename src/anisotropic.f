!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
! 
      subroutine anisotropic(anisol,anisox)
!
!     expands the 21 anisotropic elastic constants into a
!     3x3x3x3 matrix
!
      implicit none
!
      real*8 anisol(21),anisox(3,3,3,3)
!
      intent(in) anisol
!
      intent(out) anisox
!
      anisox(1,1,1,1)=anisol(1)
      anisox(1,1,1,2)=anisol(7)
      anisox(1,1,1,3)=anisol(11)
      anisox(1,1,2,1)=anisol(7)
      anisox(1,1,2,2)=anisol(2)
      anisox(1,1,2,3)=anisol(16)
      anisox(1,1,3,1)=anisol(11)
      anisox(1,1,3,2)=anisol(16)
      anisox(1,1,3,3)=anisol(4)
      anisox(1,2,1,1)=anisol(7)
      anisox(1,2,1,2)=anisol(10)
      anisox(1,2,1,3)=anisol(14)
      anisox(1,2,2,1)=anisol(10)
      anisox(1,2,2,2)=anisol(8)
      anisox(1,2,2,3)=anisol(19)
      anisox(1,2,3,1)=anisol(14)
      anisox(1,2,3,2)=anisol(19)
      anisox(1,2,3,3)=anisol(9)
      anisox(1,3,1,1)=anisol(11)
      anisox(1,3,1,2)=anisol(14)
      anisox(1,3,1,3)=anisol(15)
      anisox(1,3,2,1)=anisol(14)
      anisox(1,3,2,2)=anisol(12)
      anisox(1,3,2,3)=anisol(20)
      anisox(1,3,3,1)=anisol(15)
      anisox(1,3,3,2)=anisol(20)
      anisox(1,3,3,3)=anisol(13)
      anisox(2,1,1,1)=anisol(7)
      anisox(2,1,1,2)=anisol(10)
      anisox(2,1,1,3)=anisol(14)
      anisox(2,1,2,1)=anisol(10)
      anisox(2,1,2,2)=anisol(8)
      anisox(2,1,2,3)=anisol(19)
      anisox(2,1,3,1)=anisol(14)
      anisox(2,1,3,2)=anisol(19)
      anisox(2,1,3,3)=anisol(9)
      anisox(2,2,1,1)=anisol(2)
      anisox(2,2,1,2)=anisol(8)
      anisox(2,2,1,3)=anisol(12)
      anisox(2,2,2,1)=anisol(8)
      anisox(2,2,2,2)=anisol(3)
      anisox(2,2,2,3)=anisol(17)
      anisox(2,2,3,1)=anisol(12)
      anisox(2,2,3,2)=anisol(17)
      anisox(2,2,3,3)=anisol(5)
      anisox(2,3,1,1)=anisol(16)
      anisox(2,3,1,2)=anisol(19)
      anisox(2,3,1,3)=anisol(20)
      anisox(2,3,2,1)=anisol(19)
      anisox(2,3,2,2)=anisol(17)
      anisox(2,3,2,3)=anisol(21)
      anisox(2,3,3,1)=anisol(20)
      anisox(2,3,3,2)=anisol(21)
      anisox(2,3,3,3)=anisol(18)
      anisox(3,1,1,1)=anisol(11)
      anisox(3,1,1,2)=anisol(14)
      anisox(3,1,1,3)=anisol(15)
      anisox(3,1,2,1)=anisol(14)
      anisox(3,1,2,2)=anisol(12)
      anisox(3,1,2,3)=anisol(20)
      anisox(3,1,3,1)=anisol(15)
      anisox(3,1,3,2)=anisol(20)
      anisox(3,1,3,3)=anisol(13)
      anisox(3,2,1,1)=anisol(16)
      anisox(3,2,1,2)=anisol(19)
      anisox(3,2,1,3)=anisol(20)
      anisox(3,2,2,1)=anisol(19)
      anisox(3,2,2,2)=anisol(17)
      anisox(3,2,2,3)=anisol(21)
      anisox(3,2,3,1)=anisol(20)
      anisox(3,2,3,2)=anisol(21)
      anisox(3,2,3,3)=anisol(18)
      anisox(3,3,1,1)=anisol(4)
      anisox(3,3,1,2)=anisol(9)
      anisox(3,3,1,3)=anisol(13)
      anisox(3,3,2,1)=anisol(9)
      anisox(3,3,2,2)=anisol(5)
      anisox(3,3,2,3)=anisol(18)
      anisox(3,3,3,1)=anisol(13)
      anisox(3,3,3,2)=anisol(18)
      anisox(3,3,3,3)=anisol(6)
!
      return
      end

