!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calch0interface(nmpc,ipompc,nodempc,coefmpc,h0)
!
      implicit none
!
      integer i,j,nmpc,ist,ipompc(*),ndir,nodempc(3,*),node,index
!
      real*8 h0(3,*),coefmpc(*),h0l(3)
!
      do i=1,nmpc
         ist=ipompc(i)
         if(ist.gt.0) then
            ndir=nodempc(2,ist)
!
!           looking for MPC's tying phi between the A or A-V
!           domains and the phi-domain
!
            if(ndir.eq.5) then
               node=nodempc(1,ist)
               index=nodempc(3,ist)
               do j=1,3
                  h0l(j)=0.d0
               enddo
               if(index.ne.0) then
                  do
                     do j=1,3
                        h0l(j)=h0l(j)-coefmpc(index)*
     &                       h0(j,nodempc(1,index))
                     enddo
                     index=nodempc(3,index)
                     if(index.eq.0) exit
                  enddo
               endif
               do j=1,3
                  h0(j,node)=h0l(j)/coefmpc(ist)
               enddo
            endif
         endif
      enddo
!
      return
      end
