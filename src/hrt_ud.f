!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine hrt_ud(ielfa,vel,ipnei,nef,flux,vfa,nfacea,nfaceb,
     &           xxi,xle,gradtel,neij)
!
!     determine the facial temperature using upwind difference
!
      implicit none
!
      integer ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef,nfacea,
     &  nfaceb,neij(*)
!
      real*8 vel(nef,0:7),flux(*),vfa(0:7,*),dd,xxv(3),xxi(3,*),
     &  qp(3),xle(*),gradtel(3,*)
!
!
!
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!        unless outlet
!
c         if((iel2.le.0).and.(ielfa(3,i).ge.0)) cycle
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         if(flux(indexf).ge.0.d0) then
!
!           outflow && (neighbor || outlet)
!
            vfa(0,i)=vel(iel1,0)
         elseif(iel2.gt.0) then
!
!           inflow && neighbor
!
            vfa(0,i)=vel(iel2,0)
         endif
      enddo
!            
      return
      end
