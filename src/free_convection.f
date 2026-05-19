!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
     
      subroutine free_convection(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial,iplausi)
!          
!     Free-convection-Flow
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),iflag,iaxial,
     &     ipkon(*),kon(*),mi(*),nrhcon(*),ntmat_,nshcon(*),iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),cp,r,dvi,
     &     physcon(*),co(3,*),vold(0:mi(2),*),ttime,time,
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*)
!
      intent(in) node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,ielprop,prop,iflag,v,
     &        cp,r,physcon,dvi,set,mi,
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
      elseif ((iflag.eq.1).or.(iflag.eq.2).or.(iflag.eq.3))then
!
!    User defined flow element
!     
      endif
!
      return
      end
      

