!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_cond(imat,ntmat_,t1l,cocon,ncocon,cond)
!
      implicit none
!
!     determines the thermal conductivity
!
      integer imat,ntmat_,id,ncocon(2,*),ncoconst,seven
!
      real*8 t1l,cocon(0:6,ntmat_,*),cond
!
      seven=7
!
!     calculating the conductivity coefficients
!
      ncoconst=ncocon(1,imat)
      if(ncoconst.eq.0) then
         write(*,*) '*ERROR in materialdata_cond'
         write(*,*) 
     &        '       fluid conductivity is lacking'
         call exit(201)
      elseif(ncoconst.gt.1) then
         write(*,*) '*ERROR in materialdata_cond'
         write(*,*) 
     &        '       conductivity for fluids must be isotropic'
         call exit(201)
      endif
!     
      call ident2(cocon(0,1,imat),t1l,ncocon(2,imat),seven,id)
      if(ncocon(2,imat).eq.0) then
         cond=0.d0
         continue
      elseif(ncocon(2,imat).eq.1) then
         cond=cocon(1,1,imat)
      elseif(id.eq.0) then
         cond=cocon(1,1,imat)
      elseif(id.eq.ncocon(2,imat)) then
         cond=cocon(1,id,imat)
      else
         cond=(cocon(1,id,imat)+
     &        (cocon(1,id+1,imat)-cocon(1,id,imat))*
     &        (t1l-cocon(0,id,imat))/
     &        (cocon(0,id+1,imat)-cocon(0,id,imat)))
     &        
      endif
!     
      return
      end







