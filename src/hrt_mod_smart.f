!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine hrt_mod_smart(nface,ielfa,vel,gradtel,gamma,xlet,
     &  xxn,xxj,ipnei,betam,nef,flux,vfa)
!
!     use the modified smart scheme to determine the facial
!     temperature
!
      implicit none
!
      integer nface,ielfa(4,*),i,j,indexf,ipnei(*),iel1,iel2,nef
!
      real*8 vel(nef,0:7),gradtel(3,*),xxn(3,*),xxj(3,*),vud,vcd,
     &  gamma(*),phic,xlet(*),betam,flux(*),vfa(0:7,*)
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,ipnei,vel,vfa,flux,gradtel,xxj,xlet)
c$omp& private(i,iel2,iel1,j,indexf,vcd,vud,phic)
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
         if(flux(indexf).ge.0.d0) then
!
            vcd=vel(iel1,0)-vel(iel2,0)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel1,0))) vcd=0.d0
!
            vud=2.d0*xlet(indexf)*
     &           (gradtel(1,iel1)*xxj(1,indexf)+
     &            gradtel(2,iel1)*xxj(2,indexf)+
     &            gradtel(3,iel1)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(0,i)=vel(iel1,0)
               cycle
            endif
!     
            phic=1.d0+vcd/vud
c            write(*,*) 'calcvfa1 ',i,phic
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(0,i)=vel(iel1,0)
            elseif(phic.le.1.d0/6.d0) then
               vfa(0,i)=3.d0*vel(iel1,0)-2.d0*vel(iel2,0)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(0,i)=3.d0*vel(iel1,0)/4.d0+vel(iel2,0)/4.d0+vud/8.d0
            else
               vfa(0,i)=vel(iel1,0)/3.d0+2.d0*vel(iel2,0)/3.d0
            endif
         else
!
            vcd=vel(iel2,0)-vel(iel1,0)
            if(dabs(vcd).lt.1.d-3*dabs(vel(iel2,0))) vcd=0.d0
!
            vud=-2.d0*xlet(indexf)*
     &           (gradtel(1,iel2)*xxj(1,indexf)+
     &            gradtel(2,iel2)*xxj(2,indexf)+
     &            gradtel(3,iel2)*xxj(3,indexf))
!
            if(dabs(vud).lt.1.d-20) then
!
!           upwind difference
!
               vfa(0,i)=vel(iel2,0)
               cycle
            endif
!     
            phic=1.d0+vcd/vud
c            write(*,*) 'calcvfa2 ',i,phic
!     
            if((phic.ge.1.d0).or.(phic.le.0.d0)) then
!
!              upwind difference
!
               vfa(0,i)=vel(iel2,0)
            elseif(phic.le.1.d0/6.d0) then
               vfa(0,i)=3.d0*vel(iel2,0)-2.d0*vel(iel1,0)+2.d0*vud
            elseif(phic.le.0.7d0) then
               vfa(0,i)=3.d0*vel(iel2,0)/4.d0+vel(iel1,0)/4.d0+vud/8.d0
            else
               vfa(0,i)=vel(iel2,0)/3.d0+2.d0*vel(iel1,0)/3.d0
            endif
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
