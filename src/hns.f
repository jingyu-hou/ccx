!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
     
      subroutine hns(b,theta,rho,dg,sqrts0,xflow,h1,h2)
!
!     determine the flow depth h2 downstream of a hydraulic jump, 
!     corresponding to a upstream flow depth of h1
!     
      implicit none
!      
      real*8 b,rho,dg,sqrts0,xflow,h1,h2,c2,f,df,dh2,hk,
     &  xflow2,tth,A1,yg1,A2,yg2,dA2dh2,dyg2dh2,theta
!
      call hcrit(xflow,rho,b,theta,dg,sqrts0,hk)
!
      h2=h1*(-1.d0+dsqrt(1.d0+8.d0*(hk/h1)**3))/2.d0
!
      if(dabs(theta).lt.1.d-10) return
!
!     hns for a trapezoid, non-rectangular cross section
!
      c2=rho*rho*dg*sqrts0
      xflow2=xflow*xflow
      tth=dtan(theta)
      A1=h1*(b+h1*tth)
      yg1=h1*(3.d0*b+2.d0*h1*tth)/(6.d0*(b+h1*tth))
!
!     Newton-Raphson iterations
!
      do
         A2=h2*(b+h2*tth)
         yg2=h2*(3.d0*b+2.d0*h2*tth)/(6.d0*(b+h2*tth))
         dA2dh2=b+2.d0*h2*tth
         dyg2dh2=((3.d0*b+4.d0*h2*tth)*(b+tth)
     &                    -tth*h2*(3.d0*b+2.d0*h2*tth))/
     &                    (6.d0*(b+h2*tth)**2)
         f=A2*xflow2+c2*(A1*A1*A2*yg1-A1*A2*A2*yg2)-A1*xflow2
         df=dA2dh2*xflow2+c2*(A1*A1*yg1*dA2dh2-2.d0*A1*A2*dA2dh2*yg2
     &                       -A1*A2*A2*dyg2dh2)
         dh2=f/df
         if(dabs(dh2)/h2.lt.1.d-3) exit
         h2=h2-dh2
      enddo
!
      write(*,*) 'hns ','h1= ',h1,'h2= ',h2,'hk= ',hk
!
      return
      end
      

