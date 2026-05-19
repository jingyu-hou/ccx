!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
    
      subroutine hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
!     determine the critical depth
!     
      implicit none
!      
      real*8 xflow,rho,b,dg,sqrts0,hk,theta,tth,c1,xflow2,
     &  A,dBBdh,dAdh,BB,dhk
!
      hk=((xflow/(rho*b))**2/(dg*sqrts0))**(1.d0/3.d0)
!
      if(dabs(theta).lt.1.d-10) return
!
!     critical depth for trapezoid, non-rectangular cross section
!
      tth=dtan(theta)
      c1=rho*rho*dg*sqrts0
      xflow2=xflow*xflow
!
      do
         A=hk*(b+hk*tth)
         dBBdh=2.d0*tth
         dAdh=b+hk*dBBdh
         BB=dAdh
         dhk=(xflow2*BB-c1*A**3)/(xflow2*dBBdh-3.d0*c1*A*A*dAdh)
         if(dabs(dhk)/dhk.lt.1.d-3) exit
         hk=hk-dhk
      enddo
!
      return
      end
      

