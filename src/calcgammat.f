!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine calcgammat(ielfa,vel,gradtel,gamma,xlet,
     &  xxj,ipnei,betam,nef,flux,nfacea,nfaceb)
!
!     determine gamma for the temperature:
!        upwind difference: gamma=0
!        central difference: gamma=1
!
      implicit none
!
      integer ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef,
     &  nfacea,nfaceb
!
      real*8 vel(nef,0:7),gradtel(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*)
!
!     
!
      do i=nfacea,nfaceb
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         vcd=vel(iel2,0)-vel(iel1,0)
         if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,0))) vcd=0.d0
!
         if(flux(indexf).ge.0.d0) then
            vud=2.d0*xlet(indexf)*
     &           (gradtel(1,iel1)*xxj(1,indexf)+
     &            gradtel(2,iel1)*xxj(2,indexf)+
     &            gradtel(3,iel1)*xxj(3,indexf))
         else
            vud=2.d0*xlet(indexf)*
     &           (gradtel(1,iel2)*xxj(1,indexf)+
     &            gradtel(2,iel2)*xxj(2,indexf)+
     &            gradtel(3,iel2)*xxj(3,indexf))
         endif
!
         if(dabs(vud).lt.1.d-20) then
            gamma(i)=0.d0
            cycle
         endif
!            
         phic=1.d0-vcd/vud
!
         if(phic.ge.1.d0) then
            gamma(i)=0.d0
         elseif(phic.le.0.d0) then
            gamma(i)=0.d0
         elseif(betam.le.phic) then
            gamma(i)=1.d0
         else
            gamma(i)=phic/betam
         endif
      enddo
!            
      return
      end
