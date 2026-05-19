!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_cfd(nef,vel,shcon,nshcon,ielmatf,
     &  ntmat_,mi,cvel,vfa,cocon,ncocon,physcon,cvfa,ithermal,nface,
     &  umel,umfa,ielfa,hcfa,rhcon,nrhcon)
!
!     calculation of material properties at elements centers and
!     face centers (incompressible fluids)
!
      implicit none
!
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal(*),
     &  nshcon(2,*),nface,ncocon(2,*),ielfa(4,*),nrhcon(*)
!
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),cvel(*),vfa(0:7,*),
     &  cocon(0:6,ntmat_,*),physcon(*),cvfa(*),umel(*),umfa(*),
     &  hcfa(*),rhcon(0:1,ntmat_,*)
!
!
!
!     element (cell) values
!
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmatf(1,i)
!
!        density
!
         call materialdata_rho(rhcon,nrhcon,imat,vel(i,5),t1l,ntmat_,
     &            ithermal)
!
!        heat capacity at constant volume
!        (for liquids: =heat capacity at constant pressure)
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cvel(i),
     &       physcon)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umel(i),t1l,ntmat_,
     &            ithermal)
      enddo
!
!     facial values
!
      do i=1,nface
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmatf(1,ielfa(1,i))
!
!        density
!
         call materialdata_rho(rhcon,nrhcon,imat,vfa(5,i),t1l,ntmat_,
     &            ithermal)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cvfa(i),
     &       physcon)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umfa(i),t1l,ntmat_,
     &            ithermal)
!
!        heat conduction
!
         call materialdata_cond(imat,ntmat_,t1l,cocon,ncocon,hcfa(i))
      enddo
!            
      return
      end
