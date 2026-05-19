!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
!
!   Subroutine pre_extrapolate.f
!
!      Interpolates xstate values for the new integration points 
!      at the beginning of the new increment. 
!
!   by: Jaro Hokkanen
!
!
      subroutine interpolatestate(ne,ipkon,kon,lakon,ne0,mi,xstate,
     &  pslavsurf,nstate_,xstateini,islavsurf,islavsurfold,pslavsurfold,
     &  tieset,ntie,itiefac)
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*81 tieset(3,*)
!
      integer ipkon(*),kon(*),ne,i,n,mi(*),indexc,ne0,indexcj,
     &  nstate_,kk,nopespring,iface,ifacej,ielemslave,ll,
     &  numpts,islavsurf(2,*),islavsurfold(2,*),ntie,itiefac(2,*)
!
      real*8 xstate(nstate_,mi(1),*),pslavsurf(3,*),pslavsurfold(3,*),
     &  xstateini(nstate_,mi(1),*)
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do kk=itiefac(1,i),itiefac(2,i)
            numpts=islavsurfold(2,kk+1)-islavsurfold(2,kk)
            if(numpts.gt.2) then
               call interpolateinface(kk,xstate,xstateini,numpts,
     &              nstate_,mi,islavsurf,pslavsurf,
     &              ne0,islavsurfold,pslavsurfold)
            endif
         enddo
      enddo
!     
      return
      end
      
