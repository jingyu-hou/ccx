!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine sigini(sigma,coords,ntens,ncrds,noel,npt,layer,
     &  kspt,lrebar,rebarn)
!
!     user subroutine sigini
!
!     INPUT:
!
!     coords             coordinates of the integration point
!     ntens              number of stresses to be defined
!     ncrds              number of coordinates
!     noel               element number
!     npt                integration point number
!     layer              currently not used
!     kspt               currently not used 
!     lrebar             currently not used (value: 0)
!     rebarn             currently not used
!
!     OUTPUT:
!
!     sigma(1..ntens)    residual stress values in the integration
!                        point. If ntens=6 the order of the 
!                        components is 11,22,33,12,13,23
!           
      implicit none
!
      character*80 rebarn
      integer ntens,ncrds,noel,npt,layer,kspt,lrebar
      real*8 sigma(*),coords(*)
!
      sigma(1)=-100.d0*coords(2)
      sigma(2)=-100.d0*coords(2)
      sigma(3)=-100.d0*coords(2)
      sigma(4)=0.d0
      sigma(5)=0.d0
      sigma(6)=0.d0
!
      return
      end

