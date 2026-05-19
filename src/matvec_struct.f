!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine matvec_struct(n,x,y,nelt,ia,ja,a,isym)
      use omp_lib
!
      implicit none
!
      integer ia(*),ja(*),i,j,l,n,nelt,isym,nd,na
      real*8 y(*),x(*),a(*)
!
!     number of off-diagonal terms
!
      nd=nelt-n
!
      if(isym.eq.0) then
         na=nd/2
!
!        non-symmetric
!
!        diagonal terms
!
c$omp parallel default(none)
c$omp& shared(n,x,a,y,nd,ja,ia,na)
c$omp& private(i,l,j)
c$omp do
         do i=1,n
            y(i)=a(nd+i)*x(i)
         enddo
c$omp end do
!     
!        off-diagonal terms
!     
!        number of upper triangular terms
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
c$omp atomic
               y(i)=y(i)+a(l)*x(j)
            enddo
         enddo
c$omp end do
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
               y(j)=y(j)+a(l+na)*x(i)
            enddo
         enddo
c$omp end do
c$omp end parallel
!
      else
!
!        symmetric
!     
!        diagonal terms
!
c$omp parallel default(none)
c$omp& shared(n,x,a,y,nd,ja,ia)
c$omp& private(i,l,j)
c$omp do
         do i=1,n
            y(i)=a(nd+i)*x(i)
         enddo
c$omp end do
!
!        off-diagonal terms
!     
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
c$omp atomic
               y(i)=y(i)+a(l)*x(j)
            enddo
         enddo
c$omp end do
!
c$omp do
         do j=1,n
            do l=ja(j),ja(j+1)-1
               i=ia(l)
               y(j)=y(j)+a(l)*x(i)
            enddo
         enddo
c$omp end do
c$omp end parallel
      endif
!
      return
      end
