!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine materialdata_temp2fl(flcon,nflcon,imat,fl,
     &  t1l,ntmat_,ithermal)
!
      implicit none
!
!     determines the liquid fraction of the material from temperature
!
      integer nflcon(*),imat,two,ntmat_,id,ithermal
!
      real*8 flcon(0:1,ntmat_,*),fl,t1l
!
      two=2
!
      if(ithermal.eq.0) then
         fl=flcon(1,1,imat)
      else
         call ident2(flcon(0,1,imat),t1l,nflcon(imat),two,id)
         if(nflcon(imat).eq.0) then
            continue
         elseif(nflcon(imat).eq.1) then
            fl=flcon(1,1,imat)
         elseif(id.eq.0) then
            fl=flcon(1,1,imat)
         elseif(id.eq.nflcon(imat)) then
            fl=flcon(1,id,imat)
         else
            fl=flcon(1,id,imat)+
     &           (flcon(1,id+1,imat)-flcon(1,id,imat))*
     &           (t1l-flcon(0,id,imat))/
     &           (flcon(0,id+1,imat)-flcon(0,id,imat))
         endif
      endif
!
      return
      end
!     
