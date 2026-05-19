!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine predgmres_struct(n,b,x,nelt,ia,ja,a,isym,itol,tol,
     &  itmax,iter,err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork)
!
      implicit none
!
      integer n,nelt,ia(*),ja(*),isym,itol,itmax,iter,ierr,
     &  iunit,lrgw,igwk(*),ligw,iwork(*)
!
      real*8 b(*),x(*),a(*),tol,err,sb(*),sx(*),rgwk(*),
     &  rwork(*)
!
      external matvec_struct,msolve_struct
!
      itol=0
      tol=1.e-6
      itmax=0
      iunit=0
!
      igwk(1)=10
      igwk(2)=10
      igwk(3)=0
      igwk(4)=1
      igwk(5)=10
      ligw=20
!
      call dgmres(n,b,x,nelt,ia,ja,a,isym,matvec_struct,
     &  msolve_struct,itol,tol,itmax,
     &  iter,err,ierr,iunit,sb,sx,rgwk,lrgw,igwk,ligw,rwork,iwork)
!
      return
      end
