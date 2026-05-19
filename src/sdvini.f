!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,
     &  layer,kspt)
!
!     user subroutine sdvini
!
!
!     INPUT:
!
!     coords(1..3)       global coordinates of the integration point
!     nstatv             number of internal variables (must be
!                        defined by the user with the *DEPVAR card)
!     ncrds              number of coordinates
!     noel               element number
!     npt                integration point number
!     layer              not used
!     kspt               not used
!
!     OUTPUT:
!
!     statev(1..nstatv)  initial value of the internal state
!                        variables
!       
      implicit none
!
      integer nstatv,ncrds,noel,npt,layer,kspt,i
!
      real*8 statev(nstatv),coords(ncrds)
!
!     code for retrieving the internal state variables
!
c      do i=1,13
      do i=1,nstatv
         statev(i)=1.d0
      enddo
      return
      end

