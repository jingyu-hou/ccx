!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapolate_u(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &  ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &  vold,force,ielmat,thicke,ielprop,prop,i)
!
!     extrapolates nfield values at the integration points to the 
!     nodes for user element i
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*)
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,nfield,nk,i,ndim,
     &  iorienloc,ielorien(mi(3),*),ielmat(mi(3),*),ielprop(*)
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),orab(7,*),co(3,*),prop(*),
     &  vold(0:mi(2),*),thicke(mi(3),*)
!
      if(lakon(i)(2:3).eq.'1 ') then
         call extrapolate_u1(yi,yn,ipkon,inum,kon,lakon,nfield,nk,
     &        ne,mi,ndim,orab,ielorien,co,iorienloc,cflag,
     &        vold,force,ielmat,thicke,ielprop,prop,i)
      endif
!
      return
      end
