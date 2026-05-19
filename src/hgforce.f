!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine hgforce (fn,elas,a,gs,vl,mi,konl)
!
!     hourglass control forces for 8-node solid mean strain element
!
!     Reference: Flanagan, D.P., Belytschko, T.; "Uniform  strain hexahedron
!     and quadrilateral with orthogonal Hourglass control". Int. J. Num.
!     Meth. Engg., Vol. 17, 679-706, 1981. 
!
!     author: Otto-Ernst Bernhardi
!
      implicit none
      integer i,j,k,mi(*),konl(20)
      real*8 gs(8,4),a,elas(1),ahr
      real*8 vl(0:mi(2),20), fn(0:mi(2),*)
      real*8 hglf(3,4)
!
      ahr=elas(1)*a
c     write(6,*) "force:", ahr
!
      do i=1,3
         do k=1,4    
            hglf(i,k)=0.0d0
            do j=1,8
               hglf(i,k)=hglf(i,k)+gs(j,k)*vl(i,j)
            enddo
            hglf(i,k)=hglf(i,k)*ahr
         enddo
      enddo
      do i=1,3
         do j=1,8
            do k=1,4
               fn(i,konl(j))=fn(i,konl(j))+hglf(i,k)*gs(j,k)
            enddo
         enddo
      enddo
      return
      end
