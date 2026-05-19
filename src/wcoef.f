!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine wcoef(v,vo,al,um)
!
!     computation of the coefficients of w in the derivation of the
!     second order element stiffness matrix
!
      implicit none
!
      real*8 v(3,3,3,3),vo(3,3)
!
      real*8 a2u,al,um,au,p1,p2,p3
!
      intent(in) vo,al,um
!
      intent(inout) v
!
      a2u=al+2.d0*um
      au=al+um
!
      p1=vo(1,1)+1.d0
      p2=vo(2,2)+1.d0
      p3=vo(3,3)+1.d0
!
      v(1,1,1,1)=a2u*p1*p1+um*(vo(1,2)**2+vo(1,3)**2)
      v(2,1,1,1)=au*vo(1,2)*p1
      v(3,1,1,1)=au*vo(1,3)*p1
      v(1,2,1,1)=v(2,1,1,1)
      v(2,2,1,1)=a2u*vo(1,2)**2+um*(p1*p1+vo(1,3)**2)
      v(3,2,1,1)=au*vo(1,2)*vo(1,3)
      v(1,3,1,1)=v(3,1,1,1)
      v(2,3,1,1)=v(3,2,1,1)
      v(3,3,1,1)=a2u*vo(1,3)**2+um*(p1*p1+vo(1,2)**2)
!
      v(1,1,2,1)=al*vo(2,1)*p1+
     &  um*(2.d0*vo(2,1)*p1+vo(1,2)*p2+vo(2,3)*vo(1,3))
      v(2,1,2,1)=al*p1*p2+um*vo(2,1)*vo(1,2)
      v(3,1,2,1)=al*vo(2,3)*p1+um*vo(2,1)*vo(1,3)
      v(1,2,2,1)=al*vo(2,1)*vo(1,2)+um*p1*p2
      v(2,2,2,1)=al*vo(1,2)*p2+
     &  um*(vo(2,1)*p1+2.d0*vo(1,2)*p2+vo(2,3)*vo(1,3))
      v(3,2,2,1)=al*vo(2,3)*vo(1,2)+um*vo(1,3)*p2
      v(1,3,2,1)=al*vo(2,1)*vo(1,3)+um*vo(2,3)*p1
      v(2,3,2,1)=al*vo(1,3)*p2+um*vo(2,3)*vo(1,2)
      v(3,3,2,1)=a2u*vo(2,3)*vo(1,3)+
     &  um*(vo(2,1)*p1+vo(1,2)*p2)
!
      v(1,1,3,1)=al*vo(3,1)*p1+
     &  um*(vo(1,3)*p3+2.d0*vo(3,1)*p1+vo(3,2)*vo(1,2))
      v(2,1,3,1)=al*vo(3,2)*p1+um*vo(3,1)*vo(1,2)
      v(3,1,3,1)=al*p1*p3+um*vo(3,1)*vo(1,3)
      v(1,2,3,1)=al*vo(3,1)*vo(1,2)+um*vo(3,2)*p1
      v(2,2,3,1)=a2u*vo(3,2)*vo(1,2)+
     &  um*(vo(1,3)*p3+vo(3,1)*p1)
      v(3,2,3,1)=al*vo(1,2)*p3+um*vo(3,2)*vo(1,3)
      v(1,3,3,1)=al*vo(3,1)*vo(1,3)+um*p1*p3
      v(2,3,3,1)=al*vo(3,2)*vo(1,3)+um*vo(1,2)*p3
      v(3,3,3,1)=al*vo(1,3)*p3+
     &  um*(2.d0*vo(1,3)*p3+vo(3,1)*p1+vo(3,2)*vo(1,2))
!
      v(1,1,1,2)=al*vo(2,1)*p1+
     &  um*(vo(1,2)*p2+2.d0*vo(2,1)*p1+vo(1,3)*vo(2,3))
      v(2,1,1,2)=al*vo(1,2)*vo(2,1)+um*p1*p2
      v(3,1,1,2)=al*vo(1,3)*vo(2,1)+um*vo(2,3)*p1
      v(1,2,1,2)=al*p1*p2+um*vo(1,2)*vo(2,1)
      v(2,2,1,2)=al*vo(1,2)*p2+
     &  um*(2.d0*vo(1,2)*p2+vo(2,1)*p1+vo(1,3)*vo(2,3))
      v(3,2,1,2)=al*vo(1,3)*p2+um*vo(1,2)*vo(2,3)
      v(1,3,1,2)=al*vo(2,3)*p1+um*vo(1,3)*vo(2,1)
      v(2,3,1,2)=al*vo(1,2)*vo(2,3)+um*vo(1,3)*p2
      v(3,3,1,2)=a2u*vo(1,3)*vo(2,3)+
     &  um*(vo(1,2)*p2+vo(2,1)*p1)
!
      v(1,1,2,2)=a2u*vo(2,1)**2+um*(p2*p2+vo(2,3)**2)
      v(2,1,2,2)=au*vo(2,1)*p2
      v(3,1,2,2)=au*vo(2,3)*vo(2,1)
      v(1,2,2,2)=v(2,1,2,2)
      v(2,2,2,2)=a2u*p2*p2+um*(vo(2,1)**2+vo(2,3)**2)
      v(3,2,2,2)=au*vo(2,3)*p2
      v(1,3,2,2)=v(3,1,2,2)
      v(2,3,2,2)=v(3,2,2,2)
      v(3,3,2,2)=a2u*vo(2,3)**2+um*(p2*p2+vo(2,1)**2)
!
      v(1,1,3,2)=a2u*vo(3,1)*vo(2,1)+
     &  um*(vo(3,2)*p2+vo(2,3)*p3)
      v(2,1,3,2)=al*vo(3,2)*vo(2,1)+um*vo(3,1)*p2
      v(3,1,3,2)=al*vo(2,1)*p3+um*vo(3,1)*vo(2,3)
      v(1,2,3,2)=al*vo(3,1)*p2+um*vo(3,2)*vo(2,1)
      v(2,2,3,2)=al*vo(3,2)*p2+
     &  um*(2.d0*vo(3,2)*p2+vo(2,3)*p3+vo(3,1)*vo(2,1))
      v(3,2,3,2)=al*p2*p3+um*vo(3,2)*vo(2,3)
      v(1,3,3,2)=al*vo(3,1)*vo(2,3)+um*vo(2,1)*p3
      v(2,3,3,2)=al*vo(3,2)*vo(2,3)+um*p2*p3
      v(3,3,3,2)=al*vo(2,3)*p3+
     &  um*(vo(3,2)*p2+2.d0*vo(2,3)*p3+vo(3,1)*vo(2,1))
!
      v(1,1,1,3)=al*vo(3,1)*p1+
     &  um*(vo(1,3)*p3+2.d0*vo(3,1)*p1+vo(1,2)*vo(3,2))
      v(2,1,1,3)=al*vo(1,2)*vo(3,1)+um*vo(3,2)*p1
      v(3,1,1,3)=al*vo(1,3)*vo(3,1)+um*p1*p3
      v(1,2,1,3)=al*vo(3,2)*p1+um*vo(1,2)*vo(3,1)
      v(2,2,1,3)=a2u*vo(1,2)*vo(3,2)+
     &  um*(vo(1,3)*p3+vo(3,1)*p1)
      v(3,2,1,3)=al*vo(1,3)*vo(3,2)+um*vo(1,2)*p3
      v(1,3,1,3)=al*p1*p3+um*vo(1,3)*vo(3,1)
      v(2,3,1,3)=al*vo(1,2)*p3+um*vo(1,3)*vo(3,2)
      v(3,3,1,3)=al*vo(1,3)*p3+
     &  um*(2.d0*vo(1,3)*p3+vo(3,1)*p1+vo(1,2)*vo(3,2))
!
      v(1,1,2,3)=a2u*vo(2,1)*vo(3,1)+
     &  um*(vo(2,3)*p3+vo(3,2)*p2)
      v(2,1,2,3)=al*vo(3,1)*p2+um*vo(2,1)*vo(3,2)
      v(3,1,2,3)=al*vo(2,3)*vo(3,1)+um*vo(2,1)*p3
      v(1,2,2,3)=al*vo(2,1)*vo(3,2)+um*vo(3,1)*p2
      v(2,2,2,3)=al*vo(3,2)*p2+
     &  um*(vo(2,3)*p3+2.d0*vo(3,2)*p2+vo(2,1)*vo(3,1))
      v(3,2,2,3)=al*vo(2,3)*vo(3,2)+um*p2*p3
      v(1,3,2,3)=al*vo(2,1)*p3+um*vo(2,3)*vo(3,1)
      v(2,3,2,3)=al*p2*p3+um*vo(2,3)*vo(3,2)
      v(3,3,2,3)=al*vo(2,3)*p3+
     &  um*(2.d0*vo(2,3)*p3+vo(3,2)*p2+vo(2,1)*vo(3,1))
!
      v(1,1,3,3)=a2u*vo(3,1)**2+um*(p3*p3+vo(3,2)**2)
      v(2,1,3,3)=au*vo(3,2)*vo(3,1)
      v(3,1,3,3)=au*vo(3,1)*p3
      v(1,2,3,3)=v(2,1,3,3)
      v(2,2,3,3)=a2u*vo(3,2)**2+um*(p3*p3+vo(3,1)**2)
      v(3,2,3,3)=au*vo(3,2)*p3
      v(1,3,3,3)=v(3,1,3,3)
      v(2,3,3,3)=v(3,2,3,3)
      v(3,3,3,3)=a2u*p3*p3+um*(vo(3,1)**2+vo(3,2)**2)
!
      return
      end
