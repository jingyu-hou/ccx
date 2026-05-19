!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine linscal10(scal,konl,scall,idim,shp)
!
!     calculates a linear approximation to the quadratic interpolation
!     of the temperatures in a C3D10 element. A
!     quadratic interpolation of the temperatures leads to quadratic
!     thermal stresses, which cannot be handled by the elements 
!     displacement functions (which lead to linear stresses). Thus,
!     the temperatures are approximated by a linear function.
!
      implicit none
!
      integer konl(*),idim
!
      real*8 scal(0:idim,*),scall,shp(4,*)
!
      scall=
     &    (shp(4,1)+(shp(4,5)+shp(4,7)+shp(4,8))/2.d0)*scal(0,konl(1))
     &    +(shp(4,2)+(shp(4,5)+shp(4,6)+shp(4,9))/2.d0)*scal(0,konl(2))
     &    +(shp(4,3)+(shp(4,6)+shp(4,7)+shp(4,10))/2.d0)*scal(0,konl(3))
     &    +(shp(4,4)+(shp(4,8)+shp(4,9)+shp(4,10))/2.d0)*scal(0,konl(4))
!
      return
      end
