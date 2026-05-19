!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
    
      subroutine user_network_element(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
!     
!     user network elements
!
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(*),idirf(*),iflag,ipkon(*),kon(*),
     &     iaxial,mi(*),iplausi
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(*),R,cp,physcon(*),dvi,
     &     co(3,*),vold(0:mi(2),*),ttime,time
!
      intent(in) node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,ielprop,prop,iflag,v,
     &     cp,R,physcon,dvi,set,co,vold,mi,ttime,
     &     time,iaxial
!
      intent(inout) identity,xflow,idirf,nodef,numf,f,df,iplausi
!
!     list of different user network elements
!
!     notice that the input deck is converted into upper case when
!     being read by WeICME. So even if the user has specified "p1"
!     in his input deck, at the present stage "P1" is stored.
!
      if((lakon(nelem)(3:4).eq.'P0').or.
     &   (lakon(nelem)(3:4).eq.'0 ')) then
!
!        this just contains a skeleton file
!
         call user_network_element_p0(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
      elseif((lakon(nelem)(3:4).eq.'P1').or.
     &   (lakon(nelem)(3:4).eq.'1 ')) then
!
!        this is a working example
!
         call user_network_element_p1(node1,node2,nodem,nelem,lakon,kon,
     &     ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,dvi,numf,set,co,vold,mi,ttime,
     &     time,iaxial,iplausi)
      endif
!     
      return
      end
