!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine deltri(numpts,n,x,y,list,bin,v,e,numtri)
!     
      implicit none
!     
      integer n,i,list(*),v(3,*),e(3,*),numtri,bin(*),p,numpts
!     
      real*8 xmin,xmax,ymin,ymax,dmax,c00001,fact,x(*),y(*)
!     
      parameter(c00001=1.d0)
!     
      xmin=x(list(1))
      xmax=xmin
      ymin=y(list(1))
      ymax=ymin
      do 5 i=2,n
         p=list(i)
         xmin=min(xmin,x(p))
         xmax=max(xmax,x(p))
         ymin=min(ymin,y(p))
         ymax=max(ymax,y(p))
 5    continue
      dmax=max(xmax-xmin,ymax-ymin)
      fact=c00001/dmax
      do 10 i=1,n
         p=list(i)
         x(p)=(x(p)-xmin)*fact
         y(p)=(y(p)-ymin)*fact
 10   continue
      call bsort(n,x,y,xmin,xmax,ymin,ymax,dmax,bin,list)
      call delaun(numpts,n,x,y,list,bin,v,e,numtri)
      do 30 i=1,n
         p=list(i)
         x(p)=x(p)*dmax+xmin
         y(p)=y(p)*dmax+ymin
 30   continue
      end
      
