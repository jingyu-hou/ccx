!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine networkelementpernode(iponoel,inoel,lakon,ipkon,kon,
     &       inoelsize,nflow,ieg,ne,network)
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponoel(*),inoel(2,*),ipkon(*),kon(*),i,j,k,
     &  inoelfree,nope,indexe,node,inoelsize,nflow,ieg(*),ne,
     &  network
!
!     determining the elements belonging to the nodes of
!     the elements
!
!     network<=1: simultaneous procedure
!     network>1: alternating procedure
!
      inoelfree=1
!
      if(network.gt.1) then
         do k=1,nflow
            i=ieg(k)
            indexe=ipkon(i)
            do j=1,3
               node=kon(indexe+j)
               if(node.eq.0) cycle
               inoel(1,inoelfree)=i
               inoel(2,inoelfree)=iponoel(node)
               iponoel(node)=inoelfree
               inoelfree=inoelfree+1
            enddo
         enddo
      else
         do i=1,ne
            if(lakon(i)(1:1).eq.'D') then
               indexe=ipkon(i)
               do j=1,3
                  node=kon(indexe+j)
                  if(node.eq.0) cycle
                  inoel(1,inoelfree)=i
                  inoel(2,inoelfree)=iponoel(node)
                  iponoel(node)=inoelfree
                  inoelfree=inoelfree+1
               enddo
            endif
         enddo
      endif
!
!     size of field inoel
!
      inoelsize=inoelfree-1
!
      return
      end
