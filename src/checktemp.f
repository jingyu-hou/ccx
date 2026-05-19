!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine checktemp(t0,t1,lakon,ne,ipkon,kon)
!
!     check whether for mechanical calculations with temperature
!     conditions an initial and final temperature value was assigned
!     to each node belonging to an element
!
      implicit none
!
      character*8 lakon(*)
!
      integer ne,ipkon(*),kon(*),i,j,nope,index,node
!
      real*8 t0(*),t1(*)
!
      do i=1,ne
!
         if(lakon(i)(4:4).eq.'2') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif(lakon(i)(1:2).eq.'ES') then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
         else
            cycle
         endif
!
         index=ipkon(i)
         if(index.lt.0) cycle
!
         do j=1,nope
            node=kon(index+j)
            if(dabs(t0(node)-1.2357111317d0).lt.1.d-10) then
               write(*,*) '*ERROR in checktemp: no initial temperature'
               write(*,*) '       defined in node ',node
               call exit(201)
            elseif(dabs(t1(node)-1.2357111317d0).lt.1.d-10) then
               write(*,*) '*ERROR in checktemp: no final temperature'
               write(*,*) '       defined in node ',node
               call exit(201)
            endif
         enddo
      enddo
!
      return
      end

