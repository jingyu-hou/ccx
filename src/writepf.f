!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine writepf(d,bjr,bji,freq,nev,mode,nherm)
!
!     writes the participation factors to unit 5
!
      implicit none
!
      integer j,nev,mode,nherm
      real*8 d(*),bjr(*),bji(*),freq,pi
!
      pi=4.d0*datan(1.d0)
!
      write(5,*)
      if(mode.eq.0) then
         write(5,100) freq
      else
         write(5,101) mode
      endif
!
 100  format('P A R T I C I P A T I O N   F A C T O R S   F O R',
     &'   F R E Q U E N C Y   ',e20.13,' (CYCLES/TIME)')
 101  format('P A R T I C I P A T I O N   F A C T O R S   F O R',
     &'   M O D E   ',i5)
!
      if(nherm.eq.1) then
         write(5,*)
         write(5,*) 'MODE NO    FREQUENCY               FACTOR'
         write(5,*) '          (CYCLES/TIME)      REAL        IMAGINARY'
         write(5,*)
         do j=1,nev
            write(5,'(i7,3(2x,e14.7))') j,d(j)/(2.d0*pi),bjr(j),bji(j)
         enddo
      else
         write(5,*)
         write(5,*) 
     &'MODE NO    FREQ. (REAL)   FREQ. (IMAG)             FACTOR'
         write(5,*) 
     &'          (CYCLES/TIME)    (RAD/TIME)        REAL      IMAGINARY'
         write(5,*)
         do j=1,nev
            write(5,'(i7,4(2x,e14.7))') j,d(2*j-1)/(2.d0*pi),d(2*j),
     &            bjr(j),bji(j)
         enddo
      endif
!
      return
      end

