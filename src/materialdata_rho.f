!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_rho(rhcon,nrhcon,imat,rho,
     &  t1l,ntmat_,ithermal)
!
      implicit none
!
!     determines the density of the material
!
      integer nrhcon(*),imat,two,ntmat_,id,ithermal
!
      real*8 rhcon(0:1,ntmat_,*),rho,t1l
!
      two=2
!
      if(ithermal.eq.0) then
         rho=rhcon(1,1,imat)
      else
         call ident2(rhcon(0,1,imat),t1l,nrhcon(imat),two,id)
         if(nrhcon(imat).eq.0) then
            continue
         elseif(nrhcon(imat).eq.1) then
            rho=rhcon(1,1,imat)
         elseif(id.eq.0) then
            rho=rhcon(1,1,imat)
         elseif(id.eq.nrhcon(imat)) then
            rho=rhcon(1,id,imat)
         else
            rho=rhcon(1,id,imat)+
     &           (rhcon(1,id+1,imat)-rhcon(1,id,imat))*
     &           (t1l-rhcon(0,id,imat))/
     &           (rhcon(0,id+1,imat)-rhcon(0,id,imat))
         endif
      endif
!
      return
      end
!     
