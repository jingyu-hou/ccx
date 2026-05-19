!
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine resultsv1(nk,nactdoh,v,sol,ipompc,nodempc,coefmpc,nmpc,
     &  mi)
!
!     calculates the velocity correction (STEP 1) in the nodes
!
      implicit none
!
      integer ipompc(*),nodempc(3,*),nmpc,nk,nactdoh(0:4,*),i,j,ist,
     &  node,ndir,index,mi(*)
!
      real*8 coefmpc(*),sol(*),v(0:mi(2),*),fixed_disp
!
!     extracting the 1st velocity correction from the solution (STEP 1)
!
      do i=1,nk
         do j=1,3
            if(nactdoh(j,i).gt.0) then
               v(j,i)=sol(nactdoh(j,i))
            else
               v(j,i)=0.d0
            endif
         enddo
c         write(*,*) 'sollll ',i,(v(j,i),j=1,3)
      enddo
c      write(*,*) 'sol307',v(1,307),v(2,307),v(3,307)
!     
!     inserting the mpc information
!     
c      do i=1,nmpc
c         ist=ipompc(i)
c         node=nodempc(1,ist)
c         ndir=nodempc(2,ist)
c         index=nodempc(3,ist)
c         fixed_disp=0.d0
c         if(index.ne.0) then
c            do
c               fixed_disp=fixed_disp-coefmpc(index)*
c     &              v(nodempc(2,index),nodempc(1,index))
c               index=nodempc(3,index)
c               if(index.eq.0) exit
c            enddo
c         endif
c         fixed_disp=fixed_disp/coefmpc(ist)
c         v(ndir,node)=fixed_disp
c      enddo
!
      return
      end
