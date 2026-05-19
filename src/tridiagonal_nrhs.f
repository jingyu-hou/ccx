!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine tridiagonal_nrhs(a,b,n,m,nrhs)
!
!     solves a tridiagonal system of equations A*x=b with
!     n equations. The off-diagonals are in symmetric positions
!     about the main diagonal and m entries away; the matrix does
!     not have to be symmetric
!
!     the first line contains the elements A(1,1) and A(1,1+m)
!     the last line contains the elements A(n,n-m),A(n,n)
!     (assuming full storage mode)
!
!     tridiagonal storage: 
!       left diagonal in a(1,*)
!       main diagonal in a(2,*)
!       right diagonal in a(3,*)
!
!     INPUT:
!   
!     a                  tridiagonal matrix
!     b                  right hand side
!     n                  number of rows and columns of a
!     m                  shift of the off-diagonals compared
!                        to the main diagonal
!     nrsh               number of right hand sides
!
!     OUTPUT
!
!     b                  solution
!
      implicit none
!
      integer i,k,n,m,nrhs
!
      real*8 a(3,n),b(nrhs,n),y
!
      intent(in) n,m,nrhs
!
      intent(inout) a,b
!
!     Gaussian elimination: all values below the
!     diagonal are set to zero, the diagonal values set to 1
!
!     only the unknown values a(3,*) and b(nrhs,*) are calculated
!
      do k=1,m
         a(3,k)=a(3,k)/a(2,k)
         do i=1,nrhs
            b(i,k)=b(i,k)/a(2,k)
         enddo
      enddo
!
      do k=m+1,n-m
         y=1.d0/(a(3,k-m)*a(1,k)-a(2,k))
         a(3,k)=-a(3,k)*y
         do i=1,nrhs
            b(i,k)=(b(i,k-m)*a(1,k)-b(i,k))*y
         enddo
      enddo
!
      do k=n-m+1,n
         a(2,k)=a(3,k-m)*a(1,k)-a(2,k)
         do i=1,nrhs
            b(i,k)=(b(i,k-m)*a(1,k)-b(i,k))/a(2,k)
         enddo
      enddo
!
!     back substitution
!
      do k=n-m,1,-1
         do i=1,nrhs
            b(i,k)=b(i,k)-a(3,k)*b(i,k+m)
         enddo
      enddo
!
      return
      end

