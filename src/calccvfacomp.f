!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calccvfacomp(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
     &  mi,ielfa,cvfa,physcon)
!
!     calculation of the secant heat capacity at constant volume
!     at the face centers (compressible media)
!
      implicit none
!
      integer nface,i,nshcon(2,*),imat,ntmat_,mi(*),
     &  ielmatf(mi(3),*),ielfa(4,*)
!
      real*8 t1l,vfa(0:7,*),cp,shcon(0:3,ntmat_,*),cvfa(*),physcon(*)
!     
      do i=1,nface
         t1l=vfa(0,i)
!
!        take the material of the first adjacent element
!
         imat=ielmatf(1,ielfa(1,i))
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &       physcon)
!
!        cv=cp-r
!
         cvfa(i)=cp-shcon(3,1,imat)
      enddo
!            
      return
      end
