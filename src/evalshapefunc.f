!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine evalshapefunc(xil,etl,xl2,nopes,p)
!     
      implicit none
!     
      integer i,j,nopes,iflag
!     
      real*8  xl2(3,8),xs2(3,2),shp2(7,8),p(3),xsj2(3),
     &  xil,etl
!     
      iflag=1
      do j=1,3
         p(j)=0.d0
      enddo 
      if(nopes.eq.8)then
         call shape8q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4)then
         call shape4q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6)then
         call shape6tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      endif 
      do i=1,nopes
         do j=1,3
            p(j)=p(j)+xl2(j,i)*shp2(4,i)
         enddo
      enddo 
!     
      return
      end
