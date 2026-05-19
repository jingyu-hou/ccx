!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine gapcon(ak,d,flowm,temp,predef,time,ciname,slname,
     &   msname,coords,noel,node,npred,kstep,kinc,area)
!
!     user subroutine gapcon
!
!
!     INPUT:
!
!     d(1)               separation between the surfaces
!     d(2)               pressure transmitted across the surfaces
!     flowm              not used
!     temp(1)            temperature at the slave node (node-to-face
!                        contact) or at the slave integration point
!                        (face-to-face contact)
!     temp(2)            temperature at the corresponding master
!                        position
!     predef             not used
!     time(1)            step time at the end of the increment
!     time(2)            total time at the end of the increment
!     ciname             surface interaction name
!     slname             not used
!     msname             not used
!     coords(1..3)       coordinates of the slave node (node-to-face
!                        contact) or of the slave integration point
!                        (face-to-face contact)
!     noel               element number of the contact spring element
!     node               slave node number; zero for face-to-face contact
!     npred              not used
!     kstep              step number
!     kinc               increment number
!     area               slave area corresponding to the contact spring
!                        element
!
!     OUTPUT:
!
!     ak(1)              gap conductance
!     ak(2..5)           not used
!           
      implicit none
!
      character*80 ciname,slname,msname
!
      integer noel,node,npred,kstep,kinc
!
      real*8 ak(5),d(2),flowm(2),temp(2),predef(2,*),time(*),coords(3),
     &  area
!
      intent(in) d,flowm,temp,predef,time,ciname,slname,
     &   msname,coords,noel,node,npred,kstep,kinc,area
!
      intent(inout) ak
!
!     insert code here
!
      return
      end

