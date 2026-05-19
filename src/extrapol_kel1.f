!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_kel1(ielfa,xrlfa,vfap,vel,ifabou,
     &  nef,umfa,constant,inlet,nfacea,nfaceb)
!
!     extrapolation of turbulent kinetic energy values to the faces
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2,
     &  ipointer,inlet(*)
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),xl1,
     &  umfa(*),constant
!
!
!     
      do i=nfacea,nfaceb
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face between two elements: interpolation
!
            vfap(6,i)=xl1*vel(iel1,6)+xrlfa(2,i)*vel(iel2,6)
!
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; no zero gradient
!
!           iel2=0 is not possible: if iel2=0, there are no
!           boundary conditions on the face, hence it is an
!           exit, which means zero gradient and 
!           ielfa(3,i) <= 0
!     
            ipointer=-iel2
!     
            if(ifabou(ipointer+5).gt.0) then
!     
!              wall: kinetic turbulent energy known
!     
               vfap(6,i)=0.d0
            elseif((inlet(i).eq.1).or.
     &             (ifabou(ipointer+5).lt.0)) then
!     
!              inlet or sliding conditions: kinetic turbulent energy known
!     
               vfap(6,i)=constant*umfa(i)
            else
!
!              extrapolation
!
               vfap(6,i)=xl1*vel(iel1,6)+xrlfa(3,i)*vel(ielfa(3,i),6)
            endif
         else
!
!           boundary face; zero gradient
!
            vfap(6,i)=vel(iel1,6)
         endif
      enddo
!            
      return
      end
