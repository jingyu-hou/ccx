!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_th(cocon,ncocon,imat,iorien,pgauss,orab,
     &  ntmat_,coconloc,mattyp,t1l,rhcon,nrhcon,rho,shcon,nshcon,sph,
     &  xstiff,iint,iel,istiff,mi,xstateini,nstate_,nelcon)
!
      implicit none
!
!     determines the density, the specific heat and the conductivity 
!     in an integration point with coordinates pgauss
!
      integer ncocon(2,*),imat,iorien,k,mattyp,mi(*),
     &  ntmat_,id,two,four,seven,nrhcon(*),nshcon(*),
     &  iint,iel,ncond,istiff,ncoconst,nelcon(2,*),nstate_
!
      real*8 cocon(0:6,ntmat_,*),orab(7,*),coconloc(6),t1l,
     &  pgauss(3),rhcon(0:1,ntmat_,*),
     &  shcon(0:3,ntmat_,*),rho,sph,xstiff(27,mi(1),*),
     &  xstateini(nstate_,mi(1),*)
!
      two=2
      four=4
      seven=7
!
      if(istiff.eq.1) then
!
         ncond=ncocon(1,imat)
         if((ncond.le.-100).or.(iorien.ne.0)) ncond=6
!
!        calculating the density (needed for the capacity matrix)
!
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
!     
!        calculating the specific heat (needed for the capacity matrix)
!     
         call ident2(shcon(0,1,imat),t1l,nshcon(imat),four,id)
         if(nshcon(imat).eq.0) then
            continue
         elseif(nshcon(imat).eq.1) then
            sph=shcon(1,1,imat)
         elseif(id.eq.0) then
            sph=shcon(1,1,imat)
         elseif(id.eq.nshcon(imat)) then
            sph=shcon(1,id,imat)
         else
            sph=shcon(1,id,imat)+
     &           (shcon(1,id+1,imat)-shcon(1,id,imat))*
     &           (t1l-shcon(0,id,imat))/
     &           (shcon(0,id+1,imat)-shcon(0,id,imat))
         endif
!     
!        determining the conductivity coefficients
!
         do k=1,6
            coconloc(k)=xstiff(21+k,iint,iel)
         enddo
!     
!        determining the type: isotropic, orthotropic or anisotropic
!
         if(ncond.le.1) then
            mattyp=1
         elseif(ncond.le.3) then
            mattyp=2
         else
            mattyp=3
         endif
!
      else
!
         ncoconst=ncocon(1,imat)
         if(ncoconst.le.-100) ncoconst=-ncoconst-100
!     
!     calculating the conductivity coefficients
!
         call ident2(cocon(0,1,imat),t1l,ncocon(2,imat),seven,id)
         if(ncocon(2,imat).eq.0) then
            do k=1,6
               coconloc(k)=0.d0
            enddo
            continue
         elseif(ncocon(2,imat).eq.1) then
            do k=1,ncoconst
               coconloc(k)=cocon(k,1,imat)
            enddo
         elseif(id.eq.0) then
            do k=1,ncoconst
               coconloc(k)=cocon(k,1,imat)
            enddo
         elseif(id.eq.ncocon(2,imat)) then
            do k=1,ncoconst
               coconloc(k)=cocon(k,id,imat)
            enddo
         else
            do k=1,ncoconst
               coconloc(k)=(cocon(k,id,imat)+
     &              (cocon(k,id+1,imat)-cocon(k,id,imat))*
     &              (t1l-cocon(0,id,imat))/
     &              (cocon(0,id+1,imat)-cocon(0,id,imat)))
     &              
            enddo
         endif
      
         if (nelcon(1,imat).eq.-70) then
            coconloc(1)=coconloc(1)*(2.d0*xstateini(1,iint,iel)/
     &         (3.d0-xstateini(1,iint,iel)))
         endif
      
      endif
!
      return
      end



