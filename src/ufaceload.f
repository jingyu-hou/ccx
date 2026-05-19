!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine ufaceload(co,ipkon,kon,lakon,nboun,nodeboun,
     &  nelemload,sideload,nload,ne,nk)
!
!
!     INPUT:
!
!     co(0..3,1..nk)     coordinates of the nodes
!     ipkon(*)           element topology pointer into field kon
!     kon(*)             topology vector of all elements
!     lakon(*)           vector with elements labels
!     nboun              number of SPC's
!     nodeboun(*)        SPC node
!     nelemload(1..2,*)  1: elements faces of which are loaded
!                        2: nodes for environmental temperatures
!     sideload(*)        load label
!     nload              number of facial distributed loads
!     ne                 highest element number
!     nk                 highest node number
!
!     user routine called at the start of each step; possible use:
!     calculation of the area of sets of elements for
!     further use to calculate film or radiation coefficients.
!     The areas can be shared using common blocks.
!
      implicit none
!
      character*8 lakon(*)
      character*20 sideload(*)
!
      integer nelemload(2,*),nload,kon(*),ipkon(*),nk,ne,nboun,
     &  nodeboun(*)
!
      real*8 co(3,*)
!
!     enter code here
!
      return
      end
      


