!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_dvi(shcon,nshcon,imat,dvi,t1l,ntmat_,
     &  ithermal)
!
      implicit none
!
!     determines the dynamic viscosity
!
      integer imat,ntmat_,id,nshcon(*),four,ithermal
!
      real*8 t1l,shcon(0:3,ntmat_,*),dvi
!
      four=4
!     
      if(ithermal.eq.0) then
         dvi=shcon(2,1,imat)
      else
         call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
         if(nshcon(imat).eq.0) then
            continue
         elseif(nshcon(imat).eq.1) then
            dvi=shcon(2,1,imat)
         elseif(id.eq.0) then
            dvi=shcon(2,1,imat)
         elseif(id.eq.nshcon(imat)) then
            dvi=shcon(2,id,imat)
         else
            dvi=shcon(2,id,imat)+
     &           (shcon(2,id+1,imat)-shcon(2,id,imat))*
     &           (t1l-shcon(0,id,imat))/
     &           (shcon(0,id+1,imat)-shcon(0,id,imat))
         endif
      endif
!
      return
      end







