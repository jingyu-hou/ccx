!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine networkinum(ipkon,inum,kon,lakon,
     &  ne,itg,ntg,network)
!
!     assigns a negative sign to the inum values corresponding to
!     network nodes
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),ne,indexe,i,node1,node2,node3,
     &  itg(*),ntg,network
!
!     changing the sign of inum for all network nodes
!
      do i=1,ntg
         if(inum(itg(i)).gt.0) inum(itg(i))=-inum(itg(i))
      enddo
!
!     changing the sign of inum for all nodes belonging to 
!     "D "-elements (network elements without type in purely
!     thermal networks)
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         lakonl=lakon(i)
         if((lakonl(1:1).ne.'D ').and.
     &      ((lakonl(1:1).ne.'D').or.(network.ne.1))) cycle
!
         indexe=ipkon(i)
         if(kon(indexe+1).ne.0)  then
            node1=kon(indexe+1)
            if(inum(node1).gt.0) inum(node1)=-inum(node1)
         endif
         node2=kon(indexe+2)
         if(inum(node2).gt.0) inum(node2)=-inum(node2)
         if(kon(indexe+3).ne.0) then
            node3=kon(indexe+3)
            if(inum(node3).gt.0) inum(node3)=-inum(node3)
         endif
      enddo
!
      return
      end
