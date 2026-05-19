!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calcttel(nef,vel,shcon,nshcon,ielmatf,
     &  ntmat_,mi,physcon,ttel)
!
!     calculation of material properties at elements centers and
!     face centers (compressible fluids)
!
      implicit none
!
      integer nef,i,imat,ntmat_,mi(*),ielmatf(mi(3),*),
     &  nshcon(2,*)
!
      real*8 t1l,vel(nef,0:7),shcon(0:3,ntmat_,*),
     &  cp,physcon(*),
     &  ttel(*)
!
      intent(in) nef,shcon,nshcon,ielmatf,ntmat_,mi,physcon
!
      intent(inout) vel,ttel
!     
c$omp parallel default(none)
c$omp& shared(nef,vel,ielmatf,ntmat_,shcon,nshcon,physcon,ttel)
c$omp& private(i,t1l,imat,cp)
!
!     element (cell) values
!
c$omp do
      do i=1,nef
         t1l=vel(i,0)
         imat=ielmatf(1,i)
!     
!     heat capacity at constant volume
!     
         call materialdata_cp_sec(imat,ntmat_,t1l,shcon,nshcon,cp,
     &        physcon)
!     
         ttel(i)=t1l+(vel(i,1)*vel(i,1)+
     &        vel(i,2)*vel(i,2)+
     &        vel(i,3)*vel(i,3))/(2.d0*cp)
      enddo
c$omp end do
c$omp end parallel
!     
      return
      end
      
