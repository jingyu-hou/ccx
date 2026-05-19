!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine hgstiffness(s,elas,a,gs)
!
!     hourglass control stiffness for 8-node solid mean strain element
!
!     Reference: Flanagan, D.P., Belytschko, T.; "Uniform  strain hexahedron
!     and quadrilateral with orthogonal Hourglass control". Int. J. Num.
!     Meth. Engg., Vol. 17, 679-706, 1981. 
!
!     author: Otto-Ernst Bernhardi
!
      implicit none
!
      integer ii1,jj1,ii,jj,m1
!
      real*8 s(60,60),gs(8,4),a,elas(1),hgls,ahr
!
      intent(in) elas,a,gs
!
      intent(inout) s
!
      ahr=elas(1)*a
c     write(6,*) "stiffness:", ahr
!
      jj1=1
      do jj=1,8
         ii1=1
         do ii=1,jj
           hgls=0.0d0
           do m1=1,4
              hgls=hgls+gs(jj,m1)*gs(ii,m1)
           enddo
           hgls=hgls*ahr
           s(ii1,jj1)=s(ii1,jj1)+hgls
           s(ii1+1,jj1+1)=s(ii1+1,jj1+1)+hgls
           s(ii1+2,jj1+2)=s(ii1+2,jj1+2)+hgls
           ii1=ii1+3
         enddo
         jj1=jj1+3
      enddo
      return
      end
