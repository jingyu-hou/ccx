!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine umpc_user(x,u,f,a,jdof,n,force,iit,idiscon)
!
!     updates the coefficients in a user mpc
!
!     INPUT:
!
!     x(3,n)             Carthesian coordinates of the nodes in the
!                        user mpc.
!     u(3,n)             Actual displacements of the nodes in the
!                        user mpc.     
!     jdof               Actual degrees of freedom of the mpc terms
!     n                  number of terms in the user mpc
!     force              Actual value of the mpc force
!     iit                iteration number
!
!     OUTPUT:
!
!     f                  Actual value of the mpc. If the mpc is
!                        exactly satisfied, this value is zero
!     a(n)               coefficients of the linearized mpc
!     jdof               Corrected degrees of freedom of the mpc terms
!     idiscon            0: no discontinuity
!                        1: discontinuity
!                        If a discontinuity arises the previous
!                        results are not extrapolated at the start of
!                        a new increment
!
      implicit none
!
      integer jdof(*),n,iit,idiscon
!
      real*8 x(3,*),u(3,*),f,a(*),force
!
!
      return
      end













