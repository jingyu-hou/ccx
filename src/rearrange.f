!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine rearrange(au,irow,icol,ndim,neq)
!
!     modifies the sparse storage mode for the iterative solver of
!     Ernst Rank (pcgsolver)
!
      implicit none
!
      integer irow(*),icol(*),ndim,i,neq,k,icr,istart,idiag,kflag
      real*8 au(*)
!
      kflag=2
!
      call isortiid(irow,icol,au,ndim,kflag)
!
      istart=1
      k=irow(1)
      icr=0
      idiag=0
!
      do i=1,ndim
         if(irow(i).eq.k) then
            icr=icr+1
            cycle
         else
            call isortid(icol(istart),au(istart),icr,kflag)
            icr=1
            istart=i
            k=irow(i)
            idiag=idiag+1
            irow(idiag)=i-1
         endif
      enddo
!
!     last row
!
      call isortid(icol(istart),au(istart),icr,kflag)
      idiag=idiag+1
      irow(idiag)=ndim
!
      return
      end
