!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calcttfaext(nfaext,vfa,shcon,nshcon,ielmatf,
     &  ntmat_,mi,physcon,ttfa,ifaext,ielfa)
!
!     calculation of the total temperature at the external faces
!     based on the primary variables:
!
!     Tt=T+v**2/(2*cp)
!
      implicit none
!
      integer nfaext,i,j,imat,ntmat_,mi(*),ielmatf(mi(3),*),
     &  nshcon(2,*),ifaext(*),ielfa(4,*)
!
      real*8 t1l,vfa(0:7,*),shcon(0:3,ntmat_,*),
     &  cp,physcon(*),ttfa(*)
!
      intent(in) nfaext,shcon,nshcon,ielmatf,ntmat_,mi,physcon,ielfa
!
      intent(inout) vfa,ttfa
!     
c$omp parallel default(none)
c$omp& shared(nfaext,vfa,ielmatf,ntmat_,shcon,nshcon,physcon,ttfa,ielfa,
c$omp&        ifaext)
c$omp& private(i,j,t1l,imat,cp)
!
!     external face values
!
c$omp do
      do j=1,nfaext
         i=ifaext(j)
         t1l=vfa(0,i)
         imat=ielmatf(1,ielfa(1,i))
!     
!     heat capacity at constant volume
!     
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &        physcon)
!     
         ttfa(i)=t1l+(vfa(1,i)*vfa(1,i)+
     &        vfa(2,i)*vfa(2,i)+
     &        vfa(3,i)*vfa(3,i))/(2.d0*cp)
      enddo
c$omp end do
c$omp end parallel
!     
      return
      end
      
