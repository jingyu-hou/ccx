!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
   
!     Subroutine plane_eq.f
!
!     Creates the plane equation from three known points. 
!     Gives the z-coordinate of the fourth point as an output.
!
!     x1,y1,z1: The coordinates of the first point
!     x2,y2,z2: The coordinates of the second point
!     x3,y3,z3: The coordinates of the third point
!     x0,y0: The x and y-coordinates for the fourth point
!     output: The z-coordinate according to the x0 and y0
!
!     by: Jaro Hokkanen
!
      subroutine plane_eq(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,output)
!
      implicit none
!
      real*8 x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,output,
     &  a,b,c,d
!
      d=x1*y2*z3+y1*z2*x3+z1*x2*y3-x1*z2*y3-y1*x2*z3-z1*y2*x3
      if(d.ne.0.d0) then
         a=1.d0/d*(y2*z3+y1*z2+z1*y3-z2*y3-y1*z3-z1*y2)
      endif  
      if(d.ne.0.d0) then
         b=1.d0/d*(x1*z3+z2*x3+z1*x2-x1*z2-x2*z3-z1*x3)
      endif  
      if(d.ne.0.d0) then
         c=1.d0/d*(x1*y2+y1*x3+x2*y3-x1*y3-y1*x2-y2*x3)
      endif  
      if(d.ne.0.d0) then
         output=1.d0/c*(1.d0-a*x0-b*y0)
      else
         output=0.d0
      endif
      return
      end
