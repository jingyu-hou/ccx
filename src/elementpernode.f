!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine elementpernode(iponoel,inoel,lakon,ipkon,kon,ne)
!
      implicit none
!
      character*8 lakon(*)
!
      integer iponoel(*),inoel(2,*),ipkon(*),kon(*),i,j,ne,
     &  inoelfree,nope,indexe,node
!
      intent(in) lakon,ipkon,kon,ne
!
      intent(inout) iponoel,inoel
!
!     determining the elements belonging to the nodes of
!     the elements
!
      inoelfree=1
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).eq.'F') cycle
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif((lakon(i)(1:2).eq.'ES').or.
     &          (lakon(i)(1:2).eq.'ED')) then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
         else
            cycle
         endif
         indexe=ipkon(i)
         do j=1,nope
            node=kon(indexe+j)
            inoel(1,inoelfree)=i
            inoel(2,inoelfree)=iponoel(node)
            iponoel(node)=inoelfree
            inoelfree=inoelfree+1
         enddo
      enddo
!
      return
      end
