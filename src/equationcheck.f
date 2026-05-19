!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
     
      subroutine equationcheck(ac,nteq,nactdog,itg,ntg,nacteq,
     &  network)
!     
      implicit none
!     
      integer nactdog(0:3,*),nteq,itg(*),ntg,i,j,k,l,node,
     &  nacteq(0:3,*),network
!     
      real*8 ac(nteq,*)
!   
!     checking the columns
!
      loop1: do j=1,nteq
         do i=1,nteq
            if(dabs(ac(i,j)).gt.0.d0) cycle loop1
         enddo
         do k=1,ntg
            node=itg(k)
            do l=0,2
               if(nactdog(l,node).eq.j) then
                  if(l.eq.0) then
                     write(*,*)      '*INFO in equationcheck: temperatur
     &e in node ',node,     ' cannot be determined: probably no incoming
     & mass flow'
                  elseif(l.eq.2) then
                     write(*,*)      '*INFO in equationcheck: pressure i 
     &n node ',node,     ' cannot be determined: all incoming elements a
     &re probably critical'
                  endif
                  cycle loop1
               endif
            enddo
         enddo
      enddo loop1
!  
!     checking the rows
! 
      loop2: do i=1,nteq
         do j=1,nteq
            if(dabs(ac(i,j)).gt.0.d0) cycle loop2
         enddo
         do k=1,ntg
            node=itg(k)
            do l=0,2
               if(nacteq(l,node).eq.i) then
                  if(l.eq.0) then
                     write(*,*)      '*INFO in equationcheck: energy equ
     &ation in node ',node,     ' is identically zero: probably no incom
     &ing mass flow'
                  elseif(l.eq.2) then
                     write(*,*)      '*INFO in equationcheck: element eq 
     &uation in node',node, ' is identically zero: the element is probab
     &ly critical'
                  endif
                  cycle loop2
               endif
            enddo
         enddo
      enddo loop2
!     
      if(network.le.2) then
         write(*,*) '*ERROR in equationcheck: singular system in'
         write(*,*) '       thermal network'
         call exit(201)
      endif
!
      return
      end
      
      
