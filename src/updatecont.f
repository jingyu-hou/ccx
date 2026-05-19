!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine updatecont(koncont,ncont,co,vold,cg,straight,mi)
!
!     update geometric date of the contact master surface triangulation
!
      implicit none
!
      integer koncont(4,*),ncont,i,j,k,node,mi(*)
!
      real*8 co(3,*),vold(0:mi(2),*),cg(3,*),straight(16,*),col(3,3)
!
      do i=1,ncont
         do j=1,3
            node=koncont(j,i)
            do k=1,3
               col(k,j)=co(k,node)+vold(k,node)
            enddo
         enddo
!
!        center of gravity of the triangles
!
         do k=1,3
            cg(k,i)=col(k,1)
         enddo
         do j=2,3
            do k=1,3
               cg(k,i)=cg(k,i)+col(k,j)
            enddo
         enddo
         do k=1,3
            cg(k,i)=cg(k,i)/3.d0
         enddo
!
!        calculating the equation of the triangle plane and the planes
!        perpendicular on it and through the triangle edges
!
         call straighteq3d(col,straight(1,i))
!
      enddo
!
      return
      end
