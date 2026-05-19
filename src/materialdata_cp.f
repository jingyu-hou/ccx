!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_cp(imat,ntmat_,t1l,shcon,nshcon,cp)
!
      implicit none
!
!     determines the specific heat
!
      integer imat,ntmat_,id,nshcon(*),four
!
      real*8 t1l,shcon(0:3,ntmat_,*),cp
!
      four=4
!     
!     calculating the specific heat
!
      call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
      if(nshcon(imat).eq.0) then
         continue
      elseif(nshcon(imat).eq.1) then
         cp=shcon(1,1,imat)
      elseif(id.eq.0) then
         cp=shcon(1,1,imat)
      elseif(id.eq.nshcon(imat)) then
         cp=shcon(1,id,imat)
      else
         cp=shcon(1,id,imat)+
     &        (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &        (t1l-shcon(0,id,imat))/
     &        (shcon(0,id+1,imat)-shcon(0,id,imat))
      endif
!
      return
      end







