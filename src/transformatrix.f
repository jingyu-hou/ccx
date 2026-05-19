!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine transformatrix(xab,p,a)
!
!     determines the transformation matrix a in a point p for a carthesian 
!     (xab(7)>0) or cylindrical transformation (xab(7)<0)
!
      implicit none
!
      integer j
!
      real*8 xab(7),p(3),a(3,3),e1(3),e2(3),e3(3),dd
!
      intent(in) xab,p
!
      intent(out) a
!
      if(xab(7).gt.0) then
!
!        carthesian transformation
!
         e1(1)=xab(1)
         e1(2)=xab(2)
         e1(3)=xab(3)
!
         e2(1)=xab(4)
         e2(2)=xab(5)
         e2(3)=xab(6)
!
         dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
         do j=1,3
            e1(j)=e1(j)/dd
         enddo
!
         dd=e1(1)*e2(1)+e1(2)*e2(2)+e1(3)*e2(3)
         do j=1,3
            e2(j)=e2(j)-dd*e1(j)
         enddo
!
         dd=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
         do j=1,3
            e2(j)=e2(j)/dd
         enddo
!
         e3(1)=e1(2)*e2(3)-e2(2)*e1(3)
         e3(2)=e1(3)*e2(1)-e1(1)*e2(3)
         e3(3)=e1(1)*e2(2)-e2(1)*e1(2)
!
      else
!
!        cylindrical coordinate system in point p
!
         e1(1)=p(1)-xab(1)
         e1(2)=p(2)-xab(2)
         e1(3)=p(3)-xab(3)
!
         e3(1)=xab(4)-xab(1)
         e3(2)=xab(5)-xab(2)
         e3(3)=xab(6)-xab(3)
!
         dd=dsqrt(e3(1)*e3(1)+e3(2)*e3(2)+e3(3)*e3(3))
!
         do j=1,3
            e3(j)=e3(j)/dd
         enddo
!
         dd=e1(1)*e3(1)+e1(2)*e3(2)+e1(3)*e3(3)
!
         do j=1,3
            e1(j)=e1(j)-dd*e3(j)
         enddo
!
         dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
!
!        check whether p belongs to the cylindrical coordinate axis
!        if so, an arbitrary vector perpendicular to the axis can
!        be taken
!
         if(dd.lt.1.d-10) then
c            write(*,*) '*WARNING in transformatrix: point on axis'
            if(dabs(e3(1)).gt.1.d-10) then
               e1(2)=1.d0
               e1(3)=0.d0
               e1(1)=-e3(2)/e3(1)
            elseif(dabs(e3(2)).gt.1.d-10) then
               e1(3)=1.d0
               e1(1)=0.d0
               e1(2)=-e3(3)/e3(2)
            else
               e1(1)=1.d0
               e1(2)=0.d0
               e1(3)=-e3(1)/e3(3)
            endif
            dd=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
         endif
!
         do j=1,3
            e1(j)=e1(j)/dd
         enddo
!
         e2(1)=e3(2)*e1(3)-e1(2)*e3(3)
         e2(2)=e3(3)*e1(1)-e1(3)*e3(1)
         e2(3)=e3(1)*e1(2)-e1(1)*e3(2)
!
      endif
!
!     finding the transformation matrix
!
      do j=1,3
         a(j,1)=e1(j)
         a(j,2)=e2(j)
         a(j,3)=e3(j)
      enddo
!
      return
      end













