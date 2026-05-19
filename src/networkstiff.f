!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine networkstiff(voldl,s,imat,konl,mi,ntmat_,shcon,
     &  nshcon,rhcon,nrhcon)
!
!     calculates the stiffness of a generic networkelement
!     element label: D + blank
!
      implicit none
!
      integer konl(20),mi(*),imat,nshcon(*),nrhcon(*),ntmat_
!
      real*8 voldl(0:mi(2),9),s(60,60),gastemp,shcon(0:3,ntmat_,*),
     &  cp,r,dvi,rhcon(0:1,ntmat_,*),rho
!
      intent(in) voldl,imat,konl,mi,ntmat_,shcon,
     &  nshcon,rhcon,nrhcon
!
      intent(inout) s
!
      gastemp=(voldl(0,1)+voldl(0,3))/2.d0
!
      call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,
     &  dvi,rhcon,nrhcon,rho)
!
      if(voldl(1,2).gt.0.d0) then
         s(3,1)=-cp*voldl(1,2)
         s(3,3)=-s(3,1)
      else
         s(1,1)=-cp*voldl(1,2)
         s(1,3)=-s(1,1)
      endif
!
      return
      end

