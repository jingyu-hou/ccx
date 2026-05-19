!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writeevcscomplex(x,nx,nm,fmin,fmax)
!
!     writes the complex eigenvalues in the .dat file
!
!     nm is the nodal diameter
!
      implicit none
!
      integer j,nx,nm(nx)
      real*8 pi,fmin,fmax
      complex*16 x(nx)
!
      pi=4.d0*datan(1.d0)
!
      write(5,*)
      write(5,*) '    E I G E N V A L U E   O U T P U T'
      write(5,*)
      write(5,*) ' NODAL   MODE NO                           FREQUENCY'
      write(5,*) 'DIAMETER                     REAL PART
     &  IMAGINARY PART'
      write(5,*) '                    (RAD/TIME)      (CYCLES/TIME)
     &   (RAD/TIME)'
      write(5,*)
      do j=1,nx
         write(5,'(i5,4x,i7,3(2x,e14.7))') nm(j),j,dreal(x(j)),
     &      dreal(x(j))/(2.d0*pi),dimag(x(j))
      enddo
!
      return
      end

