!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
    
      subroutine scavenge_pump(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,ntmat_,mi,
     &        ttime,time,iaxial,iplausi)
!     
!     scavenge pump element
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),numf,node1,node2,nodem,
     &     ielprop(*),nodef(5),idirf(5),index,iflag,
     &     ipkon(*),kon(*),mi(*),ntmat_,iaxial,iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),kappa,cp,physcon(*)
     &     ,dvi,R,ttime,time
!
      intent(in) node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,ielprop,prop,iflag,v,cp,r,physcon,dvi,set,mi,
     &        ttime,time,iaxial
!
      intent(inout) identity,xflow,idirf,nodef,numf,f,df,iplausi
!  
      if (iflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif (iflag.eq.1) then      
!
      elseif (iflag.eq.2) then
!
      elseif (iflag.eq.3) then
!
      endif
      return
      end
