!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine resultsnoddir(nk,v,nactdof,b,ipompc,nodempc,coefmpc,
     &  nmpc,mi)
!
!     copying the dof-values into (idir,node) format
!     (for sensitivity results; SPC's correspond to zero)
!
      implicit none
!
      integer mi(*),nactdof(0:mi(2),*),ipompc(*),nodempc(3,*),nk,i,j,
     &  nmpc,ist,ndir,node,index
!
      real*8 v(0:mi(2),*),b(*),coefmpc(*),bnac,fixed_disp
!
      intent(in) nk,nactdof,b,ipompc,nodempc,coefmpc,
     &  nmpc,mi
!
      intent(inout) v
!
!     copying the dof-values into (idir,node) format
!     
      do i=1,nk
         do j=1,mi(2)
            if(nactdof(j,i).gt.0) then
               bnac=b(nactdof(j,i))
            else
               cycle
            endif
            v(j,i)=bnac
         enddo
      enddo
!     
!     treating the MPC's
!
      do i=1,nmpc
         ist=ipompc(i)
         node=nodempc(1,ist)
         ndir=nodempc(2,ist)
         index=nodempc(3,ist)
         fixed_disp=0.d0
         if(index.ne.0) then
            do
               fixed_disp=fixed_disp-coefmpc(index)*
     &              v(nodempc(2,index),nodempc(1,index))
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         endif
         fixed_disp=fixed_disp/coefmpc(ist)
         v(ndir,node)=fixed_disp
      enddo
!     
      return
      end
