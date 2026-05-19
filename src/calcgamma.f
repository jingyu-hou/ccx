!
!     WeICME - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!

      subroutine calcgamma(nface,ielfa,vel,gradvel,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef,flux)
!
!     determine gamma for the velocity:
!        upwind difference: gamma=0
!        central difference: gamma=1
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:7),gradvel(3,3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*),dvel1,dvel2
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,ipnei,vel,flux,gradvel,xlet,xxj,gamma,betam)
c$omp& private(i,j,iel1,iel2,indexf,dvel1,dvel2,vcd,vud,phic)
c$omp do
      do i=1,nface
         iel2=ielfa(2,i)
!
!        faces with only one neighbor need not be treated
!
         if(iel2.le.0) cycle
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
         dvel1=dsqrt(vel(iel1,1)**2+vel(iel1,2)**2+vel(iel1,3)**2)
         dvel2=dsqrt(vel(iel2,1)**2+vel(iel2,2)**2+vel(iel2,3)**2)
!
         vcd=dvel2-dvel1
!
         if(dabs(vcd).lt.1.d-3*dvel1) vcd=0.d0
!
         if(flux(indexf).ge.0.d0) then
!
            vud=2.d0*xlet(indexf)*(
     &       (vel(iel1,1)*gradvel(1,1,iel1)+
     &        vel(iel1,2)*gradvel(2,1,iel1)+
     &        vel(iel1,3)*gradvel(3,1,iel1))*xxj(1,indexf)+
     &       (vel(iel1,1)*gradvel(1,2,iel1)+
     &        vel(iel1,2)*gradvel(2,2,iel1)+
     &        vel(iel1,3)*gradvel(3,2,iel1))*xxj(2,indexf)+
     &       (vel(iel1,1)*gradvel(1,3,iel1)+
     &        vel(iel1,2)*gradvel(2,3,iel1)+
     &        vel(iel1,3)*gradvel(3,3,iel1))*xxj(3,indexf))
            vcd=vcd*dvel1
         else
            vud=2.d0*xlet(indexf)*(
     &       (vel(iel2,1)*gradvel(1,1,iel2)+
     &        vel(iel2,2)*gradvel(2,1,iel2)+
     &        vel(iel2,3)*gradvel(3,1,iel2))*xxj(1,indexf)+
     &       (vel(iel2,1)*gradvel(1,2,iel2)+
     &        vel(iel2,2)*gradvel(2,2,iel2)+
     &        vel(iel2,3)*gradvel(3,2,iel2))*xxj(2,indexf)+
     &       (vel(iel2,1)*gradvel(1,3,iel2)+
     &        vel(iel2,2)*gradvel(2,3,iel2)+
     &        vel(iel2,3)*gradvel(3,3,iel2))*xxj(3,indexf))
            vcd=vcd*dvel2
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
c$omp end do
c$omp end parallel
!            
      return
      end
