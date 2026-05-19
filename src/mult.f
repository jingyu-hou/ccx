!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine mult(matrix,trans,n)
!
      implicit none
!
      integer i,j,k,n
      real*8 matrix(3,3),trans(3,3),a(3,3)
!
!     3x3 matrix multiplication. If n=1 then
!        matrix=trans^T*matrix,
!     if n=2 then
!        matrix=matrix*trans.
!
      if(n.eq.1) then
         do i=1,3
            do j=1,3
               a(i,j)=0.d0
               do k=1,3
                  a(i,j)=a(i,j)+trans(k,i)*matrix(k,j)
               enddo
            enddo
         enddo
      elseif(n.eq.2) then
         do i=1,3
            do j=1,3
               a(i,j)=0.d0
               do k=1,3
                  a(i,j)=a(i,j)+matrix(i,k)*trans(k,j)
               enddo
            enddo
         enddo
      endif
!
      do i=1,3
         do j=1,3
            matrix(i,j)=a(i,j)
         enddo
      enddo
!
      return
      end
         
