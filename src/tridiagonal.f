!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine tridiagonal(a,b,n,m,lda)
!
!     solves a tridiagonal system of equations a*x=b with
!     n equations. The off-diagonals are in symmetric positions
!     about the main diagonal and m entries away; the matrix does
!     not have to be symmetric
!
!     the first line contains the elements a(1,1) and a(1,1+m)
!     the last line contains the elements a(n,n-m),a(n,n)
!
!     storage: 
!       left diagonal in a(*,1)
!       main diagonal in a(*,2)
!       right diagonal in a(*,3)
!
!     INPUT:
!   
!     a                  tridiagonal matrix
!     b                  right hand side
!     n                  number of rows and columns of a
!     m                  shift of the off-diagonals compared
!                        to the main diagonal
!     lda                leading dimension of a
!
!     OUTPUT
!
!     b                  solution
!
      implicit none
!
      integer k,n,m,lda
!
      real*8 a(lda,3),b(n),y
!
      intent(in) n,m,lda
!
      intent(inout) a,b
!
!     Gaussian elimination: all values below the
!     diagonal are set to zero, the diagonal values set to 1
!
!     only the unknown values a(*,3) and b(*) are calculated
!
      do k=1,m
         a(k,3)=a(k,3)/a(k,2)
         b(k)=b(k)/a(k,2)
      enddo
!
      do k=m+1,n-m
         y=1.d0/(a(k-m,3)*a(k,1)-a(k,2))
         a(k,3)=-a(k,3)*y
         b(k)=(b(k-m)*a(k,1)-b(k))*y
      enddo
!
      do k=n-m+1,n
         a(k,2)=a(k-m,3)*a(k,1)-a(k,2)
         b(k)=(b(k-m)*a(k,1)-b(k))/a(k,2)
      enddo
!
!     back substitution
!
      do k=n-m,1,-1
         b(k)=b(k)-a(k,3)*b(k+m)
      enddo
!
      return
      end

