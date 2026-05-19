!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!     
      subroutine normmpc(nmpc,ipompc,nodempc,coefmpc,inomat)
!     
!     normalizing the coefficients of MPC's for fluid applications
!     (CFD)
!     
      implicit none
!     
      integer i,node,index,nmpc,inomat(*),nodempc(3,*),ipompc(*)
!     
      real*8 coefmpc(*),size
!     
!     normalizing the MPC-coefficients
!     
      do i=1,nmpc
        index=ipompc(i)
!     
!     check whether fluid node
!     
        node=nodempc(1,index)
        if(inomat(node).eq.0) cycle
!     
!     calculating sum of the square of the MPC coefficients
!     
        size=coefmpc(index)**2
        do
          index=nodempc(3,index)
          if(index.eq.0) exit
          size=size+coefmpc(index)**2
        enddo
!
        size=dsqrt(size)
!     
!     normalizing all terms of the MPC
!     
        index=ipompc(i)
        do
          coefmpc(index)=coefmpc(index)/size
c          write(*,*) 'normmpc',i,coefmpc(index)
          index=nodempc(3,index)
          if(index.eq.0) exit
        enddo
      enddo
!     
      return
      end
      
