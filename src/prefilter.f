!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine prefilter(co,nodedesi,ndesi,xo,yo,zo,x,y,z,
     &   nx,ny,nz)               
!
!     Filtering of sensitivities      
!
      implicit none
!
      integer nodedesi(*),ndesi,m,nx(ndesi),
     &        ny(ndesi),nz(ndesi),kflag
!
      real*8 xo(ndesi),yo(ndesi),zo(ndesi),
     &       x(ndesi),y(ndesi),z(ndesi),co(3,*)
!   
!     Create set of designnodes and perform the sorting
!     needed for near3d_se
!
      do m=1,ndesi
         xo(m)=co(1,nodedesi(m))
         x(m)=xo(m)
         nx(m)=m
         yo(m)=co(2,nodedesi(m))
         y(m)=yo(m)
         ny(m)=m
         zo(m)=co(3,nodedesi(m))
         z(m)=zo(m)
         nz(m)=m
      enddo
      kflag=2
      call dsort(x,nx,ndesi,kflag)
      call dsort(y,ny,ndesi,kflag)
      call dsort(z,nz,ndesi,kflag)
!
      return        
      end
