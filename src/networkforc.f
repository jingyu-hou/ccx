!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine networkforc(vl,tnl,imat,konl,mi,ntmat_,shcon,
     &  nshcon,rhcon,nrhcon)
!
!     calculates the concentrated flux of a generic networkelement
!     element label: D + blank
!
      implicit none
!
      integer konl(20),mi(*),imat,nshcon(*),nrhcon(*),ntmat_
!
      real*8 vl(0:mi(2),20),tnl(9),gastemp,shcon(0:3,ntmat_,*),
     &  cp,r,dvi,rhcon(0:1,ntmat_,*),rho
!
      gastemp=(vl(0,1)+vl(0,3))/2.d0
!
      call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,
     &  dvi,rhcon,nrhcon,rho)
!
!     internal force = - external force
!
      if(vl(1,2).gt.0.d0) then
         tnl(3)=cp*(vl(0,3)-vl(0,1))*vl(1,2)
      else
         tnl(1)=-cp*(vl(0,1)-vl(0,3))*vl(1,2)
      endif
!
      return
      end

