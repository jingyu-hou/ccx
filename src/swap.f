!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
!     S.W. Sloan, Adv.Eng.Software,1987,9(1),34-55.
!     Permission for use with the GPL license granted by Prof. Scott
!     Sloan on 17. Nov. 2013
!
      function swap(x1,y1,x2,y2,x3,y3,xp,yp)
!     
      implicit none
!     
      real*8 x1,y1,x2,y2,x3,y3,xp,yp,x13,y13,x23,y23,
     &     x1p,y1p,x2p,y2p,cosa,cosb,sina,sinb,c00000
!     
      logical swap
!     
      parameter(c00000=0.d0)
!     
      x13=x1-x3
      y13=y1-y3
      x23=x2-x3
      y23=y2-y3
      x1p=x1-xp
      y1p=y1-yp
      x2p=x2-xp
      y2p=y2-yp
      cosa=x13*x23+y13*y23
      cosb=x2p*x1p+y1p*y2p
      if((cosa.ge.c00000).and.(cosb.ge.c00000)) then
         swap=.false.
      elseif((cosa.lt.c00000).and.(cosb.lt.c00000)) then
         swap=.true.
      else
         sina=x13*y23-x23*y13
         sinb=x2p*y1p-x1p*y2p
         if((sina*cosb+sinb*cosa).lt.c00000) then
            swap=.true.
         else
            swap=.false.
         end if
      end if
      end
