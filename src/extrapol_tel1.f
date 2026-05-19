!     
!     WeICME (Wedge Integrated Computational Materials Engineering)
!                 - A 3-dimensional finite element program.
!     
!     Developed and maintained by Shenzhen Wedge Central 
!     South Research Institute co., Ltd., Shenzhen, China
!     
!     Copy Right 2019-2023.
!
      subroutine extrapol_tel1(ielfa,xrlfa,vfap,vel,ifabou,xbounact,
     &  nef,nfacea,nfaceb)
!
!     extrapolation of temperature values to the faces
!
!     for temperature calculations the external faces have
!     as boundary condition either a specified temperature or
!     a specified flux. If the user did not apply any of these,
!     a zero specified flux is implicitly assumed
!
      implicit none
!
      integer ielfa(4,*),ifabou(*),nfacea,nfaceb,nef,i,iel1,iel2
!
      real*8 xrlfa(3,*),vfap(0:7,*),vel(nef,0:7),xbounact(*),xl1
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
            vfap(0,i)=xl1*vel(iel1,0)+xrlfa(2,i)*vel(iel2,0)
!
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; no zero gradient
!           
            if(ifabou(-iel2).gt.0) then
!
!              temperature boundary condition
!
               vfap(0,i)=xbounact(ifabou(-iel2))
            else
!
!              flux boundary condition
!
               vfap(0,i)=vel(iel1,0)
            endif
         else
!
!           boundary face; zero gradient
!
            vfap(0,i)=vel(iel1,0)
         endif
      enddo
!            
      return
      end
