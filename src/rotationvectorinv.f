!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine rotationvectorinv(c,v)
!
!     calculates rotation matrix from rotation vector
!
      implicit none
!
      real*8 c(3,3),v(3),theta,ds,dc
!
      intent(in) v
!
      intent(out) c
!
      theta=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
!
      if (theta.eq.0.d0) then
            c(1,1)=1.d0
            c(1,2)=0.d0
            c(1,3)=0.d0
            c(2,1)=0.d0
            c(2,2)=1.d0
            c(2,3)=0.d0
            c(3,1)=0.d0
            c(3,2)=0.d0
            c(3,3)=1.d0
      else      
!
            dc=dcos(theta)
            ds=dsin(theta)
!
!     C-matrix from Guido Dhondt, The Finite Element
!     Method for Three-Dimensional Thermomechanical
!     Applications p 158
!     
            c(1,1)=dc+(1.d0-dc)*v(1)*v(1)/(theta*theta)
            c(1,2)=   (1.d0-dc)*v(1)*v(2)/(theta*theta)-ds*v(3)/theta
            c(1,3)=   (1.d0-dc)*v(1)*v(3)/(theta*theta)+ds*v(2)/theta
            c(2,1)=   (1.d0-dc)*v(2)*v(1)/(theta*theta)+ds*v(3)/theta
            c(2,2)=dc+(1.d0-dc)*v(2)*v(2)/(theta*theta)
            c(2,3)=   (1.d0-dc)*v(2)*v(3)/(theta*theta)-ds*v(1)/theta
            c(3,1)=   (1.d0-dc)*v(3)*v(1)/(theta*theta)-ds*v(2)/theta
            c(3,2)=   (1.d0-dc)*v(3)*v(2)/(theta*theta)+ds*v(1)/theta
            c(3,3)=dc+(1.d0-dc)*v(3)*v(3)/(theta*theta)
      endif
!
      return
!
      end
