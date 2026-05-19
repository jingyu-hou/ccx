!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_cfd_comp1(nef,vel,shcon,nshcon,ielmatf,
     &  ntmat_,mi,cvel,physcon,ithermal,umel,nefa,nefb)
!
!     calculation of material properties at element centers
!     (compressible fluids)
!
      implicit none
!
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),ithermal(*),
     &  nshcon(2,*),nefa,nefb
!
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),cvel(*),
     &  cp,physcon(*),umel(*)
!
!
!
!     element (cell) values
!
      do i=nefa,nefb
         t1l=vel(i,0)
         imat=ielmatf(1,i)
!
!        heat capacity at constant volume
!
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvel(i)=cp-shcon(3,1,imat)
!
!        dynamic viscosity
!
         call materialdata_dvi(shcon,nshcon,imat,umel(i),t1l,ntmat_,
     &            ithermal)
      enddo
!            
      return
      end
